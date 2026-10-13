#include "../src/io/reader/parser.h"
#include "../src/constraints/transformer/constraint_transformer.h"

#include <Eigen/LU>
#include "../src/model/model.h"
#include "../src/model/solid/c3d8.h"
#include <Eigen/Geometry>

#include <gtest/gtest.h>

#include <cmath>
#include <cstring>
#include <cstdint>
#include <lz4.h>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace fem;

namespace {

std::string read_file(const std::filesystem::path& path, bool binary = false) {
    std::ifstream input(path, binary ? std::ios::binary : std::ios::in);
    return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

void replace(std::string& text, const std::string& from, const std::string& to) {
    const auto position = text.find(from);
    ASSERT_NE(position, std::string::npos) << from;
    text.replace(position, from.size(), to);
}

using Rows = std::vector<std::vector<double>>;
using Frame = std::map<std::string, Rows>;
std::vector<Frame> read_frames(const std::string& text) {
    std::vector<Frame> frames;
    std::istringstream input(text);
    std::string line, field;
    while (std::getline(input, line)) {
        if (line.rfind("LC ", 0) == 0) frames.emplace_back();
        else if (line.rfind("FIELD, NAME=", 0) == 0) {
            const auto end = line.find(',', 12);
            field = line.substr(12, end - 12);
        } else if (line == "END FIELD") field.clear();
        else if (!field.empty() && !frames.empty()) {
            std::istringstream row(line);
            std::string value;
            row >> value; // semantic node identifier
            std::vector<double> values;
            while (row >> value) values.push_back(std::strtod(value.c_str(), nullptr));
            frames.back()[field].push_back(std::move(values));
        }
    }
    return frames;
}

template<class T>
T binary_value(const std::string& data, std::size_t offset) {
    T result{};
    std::memcpy(&result, data.data() + offset, sizeof(T));
    return result;
}

void check_binary(const std::string& binary, const std::vector<Frame>& frames) {
    struct Metadata { int step; std::string name; std::uint64_t rows, components; int width; };
    std::map<std::uint64_t, Metadata> metadata;
    int checked = 0;
    for (std::size_t offset = 0; offset + 32 <= binary.size();) {
        const auto type = binary.substr(offset, 4);
        const auto compression = binary_value<std::uint32_t>(binary, offset + 4);
        const auto stored = binary_value<std::uint64_t>(binary, offset + 8);
        const auto raw = binary_value<std::uint64_t>(binary, offset + 16);
        ASSERT_LE(offset + 32 + stored, binary.size());
        std::string payload = binary.substr(offset + 32, stored);
        offset += 32 + stored;
        if (compression == 1) {
            std::string decoded(raw, '\0');
            ASSERT_EQ(LZ4_decompress_safe(payload.data(), decoded.data(), stored, raw), raw);
            payload = std::move(decoded);
        }
        if (type == "MESH") {
            EXPECT_EQ(binary_value<std::uint64_t>(payload, 0), 21u);
        } else if (type == "FMET") {
            const auto length = binary_value<std::uint16_t>(payload, 16);
            metadata[binary_value<std::uint64_t>(payload, 0)] = {
                binary_value<std::int32_t>(payload, 8), payload.substr(18, length),
                binary_value<std::uint64_t>(payload, 22 + length),
                binary_value<std::uint64_t>(payload, 30 + length),
                static_cast<unsigned char>(payload[19 + length])};
        } else if (type == "FDAT") {
            const auto& meta = metadata.at(binary_value<std::uint64_t>(payload, 0));
            if (meta.name != "PTGAP" && meta.name != "PTFORC" && meta.name != "TF") continue;
            const auto& expected = frames.at(meta.step - 1).at(meta.name);
            ASSERT_EQ(meta.rows, expected.size());
            for (std::size_t row = 0; row < meta.rows; ++row) {
                ASSERT_EQ(meta.components, expected[row].size());
                for (std::size_t component = 0; component < meta.components; ++component) {
                    const auto position = 8 + (row * meta.components + component) * meta.width;
                    const double value = meta.width == 8 ? binary_value<double>(payload, position)
                        : binary_value<float>(payload, position);
                    EXPECT_NEAR(value, expected[row][component], 1e-4);
                }
            }
            ++checked;
        }
    }
    EXPECT_EQ(checked, 9);
}

class Pretension : public testing::Test {
protected:
    std::filesystem::path directory;
    std::string deck;
    io::reader::Parser parser;

    void SetUp() override {
        directory = std::filesystem::temp_directory_path() /
            (std::string("femaster_pretension_") + testing::UnitTest::GetInstance()->current_test_info()->name());
        std::filesystem::remove_all(directory);
        std::filesystem::create_directories(directory);
        deck = read_file(std::filesystem::path(__FILE__).parent_path() / "data/pretension_face_pair_loadsteps.inp");
        // With nu=0 this is an exact uniaxial bar: gap=F L/(EA).
        replace(deck, "210000.0, 0.3", "210000.0, 0.0");
    }
    void separate_faces() {
        replace(deck, "*ELEMENT, TYPE=C3D8", "17, -1, -1, 1\n18, 1, -1, 1\n19, 1, 1, 1\n20, -1, 1, 1\n*ELEMENT, TYPE=C3D8");
        replace(deck, "2, 5, 6, 7, 8, 9, 10, 11, 12", "2, 17, 18, 19, 20, 9, 10, 11, 12");
        replace(deck, "*SURFACE, SFSET=CUT_FACE\n2, 2, 1\n*PRETENSIONSECTION, NAME=PT1, SURFACE_A=CUT_FACE",
            "*SURFACE, SFSET=CUT_A\n1, 1, 2\n*SURFACE, SFSET=CUT_B\n2, 2, 1\n*PRETENSIONSECTION, NAME=PT1, SURFACE_A=CUT_A, SURFACE_B=CUT_B");
    }
    void abaqus(bool fixed = true, bool displacement = false) {
        replace(deck, "*ELEMENT, TYPE=C3D8", "900001, 99, 88, 77\n*ELEMENT, TYPE=C3D8");
        replace(deck, "*SURFACE, SFSET=CUT_FACE\n2, 2, 1", "*SURFACE, TYPE=ELEMENT, NAME=CUT_FACE\n2, S1");
        replace(deck, "*PRETENSIONSECTION, NAME=PT1, SURFACE_A=CUT_FACE",
            "*PRE-TENSION SECTION, SURFACE=CUT_FACE, NODE=900001");
        replace(deck, "*SUPPORT, SUPPORT_COLLECTOR=BC\nENDS, 0.0, 0.0, 0.0, NAN, NAN, NAN",
            "*BOUNDARY\nENDS, 1, 3, 0");
        for (const auto* name : {"FORCE_PRETENSION", "LOCK_SOLVED_GAP", "LOCK_WITH_EXTERNAL_LOAD"}) {
            replace(deck, std::string("*LOADCASE, TYPE=LINEARSTATIC, NAME=") + name,
                std::string("*STEP, NAME=") + name + "\n*STATIC\n1, 1");
            replace(deck, "*SUPPORTS\nBC\n", "");
            replace(deck, "*END\n", "*END STEP\n");
        }
        replace(deck, "*PRETENSION, SECTION=PT1, ACTION=LOAD, CONTROL=FORCE, VALUE=1000.0",
            displacement ? "*BOUNDARY\n900001, 1, 1, 0.0035714285714285714" : "*CLOAD\n900001, 1, 1000");
        for (int i=0; i<2; ++i)
            replace(deck, "*PRETENSION, SECTION=PT1, ACTION=LOCK",
                fixed ? "*BOUNDARY, FIXED\n900001, 1, 1" : "");
    }
    void TearDown() override { std::filesystem::remove_all(directory); }
    void run(bool femr = false) {
        std::ofstream(directory / "model.inp") << deck;
        parser.run((directory / "model.inp").string(), (directory / "model").string(), {true, true, femr});
    }
    void check(double expected_gap, double expected_force = 1000.0) {
        const auto frames = read_frames(read_file(directory / "model.res"));
        ASSERT_EQ(frames.size(), 3u);
        const auto& section = *parser.model()._data->pretension_sections.at(0);
        ASSERT_EQ(section.interface_pairs.size(), 4u);
        for (std::size_t step = 0; step < frames.size(); ++step) {
            const auto& gap = frames[step].at("PTGAP");
            const auto& displacement = frames[step].at("DISPLACEMENT");
            const auto& force = frames[step].at("PTFORC");
            double resultant = 0;
            for (const auto& pair : section.interface_pairs) {
                ASSERT_LT(static_cast<std::size_t>(pair.side_b), gap.size());
                EXPECT_NEAR(gap.at(pair.side_a).at(0), expected_gap, 1e-8);
                EXPECT_NEAR(gap.at(pair.side_b).at(0), expected_gap, 1e-8);
                Vec3 relative;
                for (Dim component = 0; component < 3; ++component) {
                    relative(component) = displacement.at(pair.side_b).at(component)
                        - displacement.at(pair.side_a).at(component);
                    resultant += force.at(pair.side_b).at(component) * section.axis_direction(component);
                }
                EXPECT_NEAR(relative.dot(section.axis_direction), expected_gap, 1e-8);
                EXPECT_NEAR((relative - expected_gap * section.axis_direction).norm(), 0, 1e-8);
            }
            EXPECT_NEAR(resultant, expected_force + (step == 2 ? section.axis_direction(2) * 4.0 / 3.0 : 0.0), 1e-3);
        }
        const auto dofs = parser.model().build_structural_dof_index_matrix();
        EXPECT_GE(dofs(section.control_node, 0), 0);
        for (Dim c = 1; c < 6; ++c) EXPECT_EQ(dofs(section.control_node, c), -1);
        for (const auto& frame : frames) {
            EXPECT_NEAR(frame.at("DISPLACEMENT").at(section.control_node).at(0), expected_gap, 1e-8);
            EXPECT_NEAR(frame.at("PTGAP").at(section.control_node).at(0), expected_gap, 1e-8);
        }
        EXPECT_NEAR(section.total_force, frames.back().at("PTFORC").at(section.control_node).at(0), 1e-3);
        if (!section.nodal_control) EXPECT_NEAR(section.locked_gap, expected_gap, 1e-10);
        EXPECT_NE(read_file(directory / "model.frd").find("PTGAP"), std::string::npos);
    }
};

TEST_F(Pretension, ForceLockExternalLoadAndBinaryOutput) {
    run(true);
    check(1000.0 * 3.0 / (210000.0 * 4.0));
    EXPECT_EQ(parser.model()._data->positions->rows, 21u);
    const auto binary = read_file(directory / "model.femr", true);
    check_binary(binary, read_frames(read_file(directory / "model.res")));
}

TEST_F(Pretension, DisplacementLockExternalLoad) {
    replace(deck, "CONTROL=FORCE, VALUE=1000.0", "CONTROL=DISPLACEMENT, VALUE=0.0035714285714285714");
    run();
    check(1000.0 * 3.0 / (210000.0 * 4.0));
}

TEST_F(Pretension, EliminationForceLockExternalLoad) {
    std::size_t position = 0;
    while ((position = deck.find("TYPE=NULLSPACE", position)) != std::string::npos) {
        deck.replace(position, 14, "TYPE=ELIMINATION");
        position += 16;
    }
    run();
    check(1000.0 * 3.0 / (210000.0 * 4.0));
}

TEST_F(Pretension, RejectsLagrange) {
    // Reject both an actively loaded section and an unloaded (tied) section.
    for (const bool loaded : {true, false}) {
        auto original = deck;
        replace(deck, "TYPE=NULLSPACE", "TYPE=LAGRANGE");
        if (!loaded) replace(deck,
            "*PRETENSION, SECTION=PT1, ACTION=LOAD, CONTROL=FORCE, VALUE=1000.0", "");
        try {
            run();
            FAIL() << "LAGRANGE must be rejected before solving";
        } catch (const std::exception& error) {
            EXPECT_NE(std::string(error.what()).find("Use NULLSPACE or ELIMINATION"), std::string::npos);
        }
        deck = original;
    }
}

TEST_F(Pretension, LagrangeWithoutPretensionKeepsExistingBehavior) {
    replace(deck, "*PRETENSIONSECTION, NAME=PT1, SURFACE_A=CUT_FACE", "");
    replace(deck, "*PRETENSION, SECTION=PT1, ACTION=LOAD, CONTROL=FORCE, VALUE=1000.0", "");
    for (int step = 0; step < 2; ++step)
        replace(deck, "*PRETENSION, SECTION=PT1, ACTION=LOCK", "");
    for (int step = 0; step < 3; ++step)
        replace(deck, "TYPE=NULLSPACE", "TYPE=LAGRANGE");
    run();
    EXPECT_TRUE(parser.model()._data->pretension_sections.empty());
    const auto frames = read_frames(read_file(directory / "model.res"));
    ASSERT_EQ(frames.size(), 3u);
    for (const auto& frame : frames) {
        EXPECT_EQ(frame.count("PTGAP"), 0u);
        EXPECT_EQ(frame.count("PTFORC"), 0u);
        for (const auto& row : frame.at("DISPLACEMENT"))
            for (const auto value : row) EXPECT_TRUE(std::isfinite(value));
    }
    for (const auto& row : frames[0].at("DISPLACEMENT"))
        for (const auto value : row) EXPECT_NEAR(value, 0, 1e-12);
    EXPECT_LT(frames[2].at("DISPLACEMENT").at(8).at(2), 0);
}

class PretensionExact : public Pretension, public testing::WithParamInterface<const char*> {
protected:
    void select_method() {
        std::size_t position = 0;
        const std::string method = std::string("TYPE=") + GetParam();
        while ((position = deck.find("TYPE=NULLSPACE", position)) != std::string::npos) {
            deck.replace(position, 14, method);
            position += method.size();
        }
    }
    void add_equations(const std::string& equations) {
        replace(deck, "*LOADCASE, TYPE=LINEARSTATIC, NAME=FORCE_PRETENSION",
            equations + "*LOADCASE, TYPE=LINEARSTATIC, NAME=FORCE_PRETENSION");
    }
};

TEST_P(PretensionExact, ForceAndDisplacementWithMixedRedundantConstraints) {
    const auto original = deck;
    for (const bool displacement : {false, true}) {
        deck = original;
        select_method();
        if (displacement) replace(deck, "CONTROL=FORCE, VALUE=1000.0",
            "CONTROL=DISPLACEMENT, VALUE=0.0035714285714285714");
        // A manual transverse MPC away from the cut, its scaled duplicate,
        // and an equation already implied by the end supports.
        add_equations("*EQUATION\n2\n9, 1, 1, 10, 1, -1\n"
                      "*EQUATION\n2\n9, 1, 2, 10, 1, -2\n"
                      "*EQUATION\n2\n1, 3, 1, 13, 3, -1\n");
        run();
        check(1000.0 * 3.0 / (210000.0 * 4.0));
    }
}

TEST_P(PretensionExact, RejectsContradictoryMixedConstraints) {
    select_method();
    replace(deck, "ENDS, 0.0, 0.0, 0.0, NAN, NAN, NAN",
        "ENDS, 0.0, 0.0, 0.0, NAN, NAN, NAN\n13, 0.0, 0.0, 0.01, NAN, NAN, NAN");
    add_equations("*EQUATION\n2\n1, 3, 1, 13, 3, -1\n");
    try {
        run();
        FAIL() << "Conflicting supports and equations must be rejected";
    } catch (const std::exception& error) {
        EXPECT_NE(std::string(error.what()).find("inconsistent"), std::string::npos);
    }
    EXPECT_FALSE(parser.model()._data->pretension_sections.at(0)->has_solved_gap);
}

TEST_P(PretensionExact, RejectsConstraintsConflictingWithPrescribedOrLockedGap) {
    const auto original = deck;
    for (const bool locked : {false, true}) {
        deck = original;
        select_method();
        separate_faces();
        const std::string equation = "*EQUATION\n2\n5, 3, 1, 17, 3, -1\n";
        if (locked) {
            // EQUATION definitions are global. Activate new supports only in
            // the LOCK step to conflict with the previously solved nonzero gap.
            replace(deck, "*LOADCASE, TYPE=LINEARSTATIC, NAME=LOCK_SOLVED_GAP",
                "*LOADCASE, TYPE=LINEARSTATIC, NAME=LOCK_SOLVED_GAP\n"
                "*SUPPORT\n5, NAN, NAN, 0, NAN, NAN, NAN\n"
                "17, NAN, NAN, 0, NAN, NAN, NAN");
        } else {
            replace(deck, "CONTROL=FORCE, VALUE=1000.0",
                "CONTROL=DISPLACEMENT, VALUE=0.0035714285714285714");
            add_equations(equation);
        }
        try {
            run();
            FAIL() << "Zero-gap constraints cannot coexist with a nonzero prescribed gap";
        } catch (const std::exception& error) {
            EXPECT_NE(std::string(error.what()).find("inconsistent"), std::string::npos);
        }
        const auto& section = *parser.model()._data->pretension_sections.at(0);
        EXPECT_EQ(section.has_solved_gap, locked);
        if (locked) EXPECT_NEAR(std::abs(section.last_solved_gap),
            1000.0 * 3.0 / (210000.0 * 4.0), 1e-10);
    }
}

INSTANTIATE_TEST_SUITE_P(Methods, PretensionExact,
    testing::Values("NULLSPACE", "ELIMINATION"));

TEST(PretensionConstraintDecision, LagrangeRegularizationAndRedundantRows) {
    using Transformer = constraint::ConstraintTransformer;
    SystemDofIds dofs(1, 6);
    dofs.setConstant(-1);
    dofs(0, 0) = 0;
    // Same DOF prescribed twice: admissible for reduced methods, singular
    // for exact multipliers. No pretension sections are involved here.
    constraint::Equations equations = {{{{0, 0, 1}}, 0}, {{{0, 0, 1}}, 0}};
    SparseMatrix stiffness(1, 1);
    stiffness.insert(0, 0) = 100;
    Transformer::Options options;
    options.method = Transformer::Method::Lagrange;
    Transformer regularized(equations, dofs, 1, options);
    const auto matrix = regularized.assemble_system_matrix(stiffness);
    EXPECT_NEAR(matrix.coeff(1, 1), -1e-8, 1e-20);
    EXPECT_NEAR(matrix.coeff(2, 2), -1e-8, 1e-20);
    EXPECT_FALSE(regularized.rank_known());
    options.lagrange_regularization = 0;
    Transformer exact(equations, dofs, 1, options);
    const DynamicMatrix dense(exact.assemble_system_matrix(stiffness));
    EXPECT_LT(dense.fullPivLu().rank(), dense.rows());
}

TEST_F(Pretension, SeparateMatchingFaces) {
    separate_faces();
    // SURFACE_A points in the opposite direction in this equivalent deck.
    replace(deck, "VALUE=1000.0", "VALUE=-1000.0");
    run();
    const auto& section = *parser.model()._data->pretension_sections.at(0);
    EXPECT_EQ(parser.model()._data->positions->rows, 21u);
    check(-1000.0 * 3.0 / (210000.0 * 4.0), -1000.0);
}

TEST_F(Pretension, ExplicitNativeControlNode) {
    replace(deck, "*ELEMENT, TYPE=C3D8", "900001, 99, 88, 77\n*ELEMENT, TYPE=C3D8");
    replace(deck, "NAME=PT1, SURFACE_A=CUT_FACE", "NAME=PT1, SURFACE_A=CUT_FACE, NODE=900001");
    run();
    check(1000.0 * 3.0 / (210000.0 * 4.0));
    EXPECT_EQ(parser.model()._data->pretension_sections[0]->control_node,
        parser.model().compiled_node_id("900001"));
    EXPECT_NE(read_file(directory / "model.res").find("control   : 900001 (DOF 1)"), std::string::npos);
}

TEST_F(Pretension, RejectsStructuralControlNode) {
    replace(deck, "NAME=PT1, SURFACE_A=CUT_FACE", "NAME=PT1, SURFACE_A=CUT_FACE, NODE=5");
    EXPECT_THROW(run(), std::exception);
}


TEST_P(PretensionExact, AbaqusNativeEquivalenceAndFixed) {
    const auto original = deck;
    for (const bool displacement : {false, true}) {
        deck = original;
        select_method();
        if (displacement) replace(deck, "CONTROL=FORCE, VALUE=1000.0",
            "CONTROL=DISPLACEMENT, VALUE=0.0035714285714285714");
        run();
        const auto native = read_frames(read_file(directory / "model.res"));
        deck = original;
        select_method();
        abaqus(true, displacement);
        run(true);
        check(1000.0 * 3.0 / (210000.0 * 4.0));
        const auto result = read_frames(read_file(directory / "model.res"));
        const auto& section = *parser.model()._data->pretension_sections[0];
        // Map all physical nodes; generated cut nodes shift by one because
        // the explicit Abaqus control node precedes mesh preparation.
        for (std::size_t step = 0; step < native.size(); ++step) {
            for (const auto* field : {"DISPLACEMENT", "PTFORC", "PTGAP"})
                for (std::size_t node=0; node<20; ++node)
                    for (std::size_t c=0; c<native[step].at(field)[node].size(); ++c)
                        EXPECT_NEAR(result[step].at(field)[node < 16 ? node : node+1][c],
                            native[step].at(field)[node][c], 1e-5);
            const double total = result[step].at("TF")[section.control_node][0];
            const double reaction = result[step].at("REACTION_FORCES")[section.control_node][0];
            EXPECT_NEAR(total, reaction + (displacement ? 0 : 1000), 1e-4);
            EXPECT_NEAR(total, result[step].at("PTFORC")[section.control_node][0], 1e-4);
        }
        EXPECT_NE(read_file(directory / "model.frd").find("TF"), std::string::npos);
        check_binary(read_file(directory / "model.femr", true), result);
    }
}

TEST_F(Pretension, AbaqusForceCarriesWithoutFixed) {
    abaqus(false);
    run();
    const auto frames = read_frames(read_file(directory / "model.res"));
    const ID q = parser.model()._data->pretension_sections[0]->control_node;
    for (const auto& frame : frames) EXPECT_NEAR(frame.at("TF")[q][0], 1000, 1e-4);
    EXPECT_GT(std::abs(frames[2].at("PTGAP")[q][0] - frames[0].at("PTGAP")[q][0]), 1e-7);
}

TEST_F(Pretension, AbaqusModNewAndAmplitude) {
    abaqus(false);
    replace(deck, "*PRE-TENSION SECTION", "*NSET, NSET=QALIAS\n900001\n*PRE-TENSION SECTION");
    replace(deck, "*STEP, NAME=FORCE_PRETENSION",
        "*AMPLITUDE, NAME=RAMP\n0, 0, 1, 0.5\n*STEP, NAME=FORCE_PRETENSION");
    replace(deck, "*CLOAD\n900001, 1, 1000", "*CLOAD, AMPLITUDE=RAMP\n900001, 1, 1000");
    replace(deck, "*STEP, NAME=LOCK_SOLVED_GAP\n*STATIC\n1, 1",
        "*STEP, NAME=LOCK_SOLVED_GAP\n*STATIC\n1, 1\n*CLOAD, OP=MOD\nQALIAS, 1, 2000");
    replace(deck, "*LOADS\nEXTERNAL", "*CLOAD, OP=NEW\nMIDDLE, 3, -1");
    run();
    const auto frames = read_frames(read_file(directory / "model.res"));
    const ID q = parser.model()._data->pretension_sections[0]->control_node;
    EXPECT_NEAR(frames[0].at("TF")[q][0], 500, 1e-4);
    EXPECT_NEAR(frames[1].at("TF")[q][0], 2000, 1e-4);
    EXPECT_NEAR(frames[2].at("TF")[q][0], 0, 1e-4);
}

TEST_F(Pretension, AbaqusBoundaryNewReleasesFixed) {
    abaqus();
    replace(deck, "*LOADS\nEXTERNAL", "");
    // Only the second FIXED occurrence (in step 3) is replaced.
    auto pos = deck.rfind("*BOUNDARY, FIXED\n900001, 1, 1");
    ASSERT_NE(pos, std::string::npos);
    deck.replace(pos, std::string("*BOUNDARY, FIXED\n900001, 1, 1").size(),
        "*BOUNDARY, OP=NEW\nENDS, 1, 3, 0\n*CLOAD, OP=MOD\n900001, 1, 2000");
    run();
    const auto frames = read_frames(read_file(directory / "model.res"));
    const ID q = parser.model()._data->pretension_sections[0]->control_node;
    EXPECT_NEAR(frames[2].at("TF")[q][0], 2000, 1e-4);
    EXPECT_NEAR(frames[2].at("PTGAP")[q][0], 2000.0*3/(210000*4), 1e-8);
}

TEST_F(Pretension, AbaqusSingletonAndNormal) {
    abaqus();
    replace(deck, "*PRE-TENSION SECTION, SURFACE=CUT_FACE, NODE=900001",
        "*NSET, NSET=CONTROL\n900001\n*PRE-TENSION SECTION, SURFACE=CUT_FACE, NODE=CONTROL, FOLLOWER=NO\n0, 0, -2");
    run();
    check(1000.0*3/(210000*4));
}

TEST_F(Pretension, AbaqusUnloadedSectionKeepsNativeZeroGap) {
    abaqus(false);
    replace(deck, "*CLOAD\n900001, 1, 1000", "");
    run();
    const auto frames = read_frames(read_file(directory / "model.res"));
    const ID q = parser.model()._data->pretension_sections[0]->control_node;
    for (const auto& frame : frames) EXPECT_NEAR(frame.at("PTGAP")[q][0], 0, 1e-10);
}

TEST_F(Pretension, AbaqusRejectsUnsupportedAndInvalidControls) {
    abaqus();
    const auto original = deck;
    const std::vector<std::pair<std::string,std::string>> cases = {
        {"NODE=900001", "NODE=999999"},
        {"NODE=900001", "NODE=5"},
        {"NODE=900001", "NODE=ENDS"},
        {"NODE=900001", "NODE=900001, FOLLOWER=YES"},
        {"SURFACE=CUT_FACE", "ELEMENT=2"},
        {"NODE=900001", "NODE=900001\n0, 0, 0"},
        {"900001, 1, 1000", "900001, 2, 1000"},
        {"900001, 1, 1000", "900001, 2, 0"},
        {"900001, 1, 1000", "900001, 1, NAN"},
        {"900001, 1, 1000", "900001, 1, 1000\n*BOUNDARY\n900001, 2, 2, 0"},
        {"900001, 1, 1000", "900001, 1, 1000\n*BOUNDARY\n900001, 1, 1, NAN"},
        {"*PRE-TENSION SECTION, SURFACE=CUT_FACE, NODE=900001",
         "*NSET, NSET=QALIAS\n900001\n*PRE-TENSION SECTION, SURFACE=CUT_FACE, NODE=900001\n*PRE-TENSION SECTION, SURFACE=CUT_FACE, NODE=QALIAS"},
        {"*PRE-TENSION SECTION, SURFACE=CUT_FACE, NODE=900001",
         "*PART, NAME=UNSUPPORTED\n*PRE-TENSION SECTION, SURFACE=CUT_FACE, NODE=900001\n*END PART"},
        {"*CLOAD\n900001, 1, 1000", "*BOUNDARY, FIXED\n900001, 1, 1"},
        {"*PRE-TENSION SECTION, SURFACE=CUT_FACE, NODE=900001",
         "*PRE-TENSION SECTION, SURFACE=CUT_FACE, NODE=900001\n*PRE-TENSION SECTION, SURFACE=CUT_FACE, NODE=900001"},
        {"*STEP, NAME=FORCE_PRETENSION",
         "*NSET, NSET=QNODE\n900001\n*TRANSFORM, NSET=QNODE, TYPE=R\n1,0,0,0,1,0\n*STEP, NAME=FORCE_PRETENSION"}
    };
    for (const auto& [from,to] : cases) {
        SCOPED_TRACE(to);
        deck = original;
        replace(deck, from, to);
        EXPECT_THROW(run(), std::exception);
    }
}

TEST_F(Pretension, AbaqusRejectsConflictingAdjustment) {
    abaqus();
    replace(deck, "*CLOAD\n900001, 1, 1000",
        "*BOUNDARY\n900001, 1, 1, 0.001\n*BOUNDARY\nNALL, 1, 1, 0");
    EXPECT_THROW(run(), std::exception);
}


TEST_P(PretensionExact, TwoIndependentRotatedSections) {
    const auto original = read_file(std::filesystem::path(__FILE__).parent_path() / "data/pretension_two_sections.inp");
    std::vector<Frame> native;
    for (const bool nodal : {false, true}) {
        deck = original;
        select_method();
        if (nodal) {
            abaqus();
            replace(deck, "*ELEMENT, TYPE=C3D8", "900002, 0, 0, 0\n*ELEMENT, TYPE=C3D8");
            replace(deck, "*PRETENSIONSECTION, NAME=PT2, SURFACE_A=CUT_SECOND",
                "*PRE-TENSION SECTION, SURFACE=CUT_SECOND, NODE=900002");
            replace(deck, "*PRETENSION, SECTION=PT2, ACTION=LOAD, CONTROL=DISPLACEMENT, VALUE=0.007142857142857143",
                "*BOUNDARY\n900002, 1, 1, 0.007142857142857143");
            for (int i=0; i<2; ++i) replace(deck, "*PRETENSION, SECTION=PT2, ACTION=LOCK",
                "*BOUNDARY, FIXED\n900002, 1, 1");
        }
        run();
        const auto frames = read_frames(read_file(directory / "model.res"));
        const auto& sections = parser.model()._data->pretension_sections;
        ASSERT_EQ(sections.size(), 2u);
        EXPECT_NE(sections[0]->control_node, sections[1]->control_node);
        EXPECT_NEAR(sections[1]->axis_direction(0), -1, 1e-12);
        for (const auto& frame : frames) {
            EXPECT_NEAR(frame.at("PTGAP")[sections[0]->control_node][0], 1000.0*3/(210000*4), 1e-8);
            EXPECT_NEAR(frame.at("PTGAP")[sections[1]->control_node][0], 2000.0*3/(210000*4), 1e-8);
            EXPECT_NEAR(frame.at("TF")[sections[1]->control_node][0], 2000, 1e-4);
        }
        if (!nodal) native = frames;
        else for (std::size_t step=0; step<3; ++step)
            for (std::size_t node=0; node<32; ++node)
                for (std::size_t c=0; c<3; ++c)
                    EXPECT_NEAR(frames[step].at("DISPLACEMENT")[node][c],
                        native[step].at("DISPLACEMENT")[node][c], 1e-8);
    }
}

TEST_F(Pretension, AbaqusAssemblyWithInstanceQualifiedReferences) {
    abaqus();
    // Move the topology and regions to a reusable Part; declare the section in
    // Assembly using compiled instance-qualified node and surface references.
    const std::string material = "*MATERIAL, NAME=STEEL\n*ELASTIC, TYPE=ISOTROPIC\n210000.0, 0.0\n";
    replace(deck, material, "");
    replace(deck, "*NODE, NSET=NALL", material + "*PART, NAME=BOLT\n*NODE, NSET=NALL");
    replace(deck, "*PRE-TENSION SECTION, SURFACE=CUT_FACE, NODE=900001\n", "");
    replace(deck, "*BOUNDARY\nENDS, 1, 3, 0",
        "*END PART\n*ASSEMBLY, NAME=A\n*INSTANCE, NAME=B1, PART=BOLT\n*END INSTANCE\n"
        "*PRE-TENSION SECTION, SURFACE=B1.CUT_FACE, NODE=B1.900001\n*END ASSEMBLY\n"
        "*BOUNDARY\nB1.ENDS, 1, 3, 0");
    replace(deck, "MIDDLE, 0.0", "B1.MIDDLE, 0.0");
    replace(deck, "900001, 1, 1000", "B1.900001, 1, 1000");
    for (int i=0; i<2; ++i) replace(deck, "\n900001, 1, 1\n", "\nB1.900001, 1, 1\n");
    run();
    check(1000.0*3/(210000*4));
}

TEST_F(Pretension, AbaqusReversedNormalAndBoundaryAmplitude) {
    abaqus(true, true);
    replace(deck, "NODE=900001\n", "NODE=900001\n0, 0, 1\n");
    replace(deck, "*STEP, NAME=FORCE_PRETENSION",
        "*AMPLITUDE, NAME=HALF\n0, 0, 1, 0.5\n*STEP, NAME=FORCE_PRETENSION");
    replace(deck, "*BOUNDARY\n900001, 1, 1, 0.0035714285714285714",
        "*BOUNDARY, AMPLITUDE=HALF\n900001, 1, 1, 0.007142857142857143");
    run();
    check(1000.0*3/(210000*4));
    EXPECT_NEAR(parser.model()._data->pretension_sections[0]->axis_direction(2), 1, 1e-12);
}

TEST_F(Pretension, RejectsBalancingOfScalarControlLoads) {
    const auto original = deck;
    for (const auto* command : {"INERTIARELIEF", "REBALANCELOADS"}) {
        deck = original;
        replace(deck, "*PRETENSION, SECTION=PT1, ACTION=LOAD", std::string("*") + command + "\n*PRETENSION, SECTION=PT1, ACTION=LOAD");
        EXPECT_THROW(run(), std::exception);
    }
}

TEST_F(Pretension, AbaqusBoundaryAliasReplacesScalarPrescription) {
    abaqus(true, true);
    replace(deck, "*PRE-TENSION SECTION", "*NSET, NSET=QALIAS\n900001\n*PRE-TENSION SECTION");
    replace(deck, "*BOUNDARY, FIXED\n900001, 1, 1",
        "*BOUNDARY, OP=MOD\nQALIAS, 1, 1, 0.007142857142857143");
    run();
    const auto frames = read_frames(read_file(directory / "model.res"));
    const ID q = parser.model()._data->pretension_sections[0]->control_node;
    EXPECT_NEAR(frames[0].at("TF")[q][0], 1000, 1e-4);
    EXPECT_NEAR(frames[1].at("TF")[q][0], 2000, 1e-4);
    EXPECT_NEAR(frames[2].at("PTGAP")[q][0], 2000.0*3/(210000*4), 1e-8);
}

TEST_F(Pretension, RejectsLockBeforeLoad) {
    replace(deck, "ACTION=LOAD, CONTROL=FORCE, VALUE=1000.0", "ACTION=LOCK");
    EXPECT_THROW(run(), std::exception);
}

TEST_F(Pretension, RejectsExteriorFaceAsMergedCut) {
    replace(deck, "2, 2, 1\n*PRETENSIONSECTION", "2, 1, 1\n*PRETENSIONSECTION");
    EXPECT_THROW(run(), std::exception);
    EXPECT_EQ(parser.model()._data->positions->rows, 16);
    EXPECT_EQ(parser.model()._data->node_mapping.size(), 16);
    EXPECT_EQ(parser.model()._data->pretension_sections[0]->control_node, -1);
}

TEST_F(Pretension, RejectsMissingLoadValue) {
    replace(deck, ", VALUE=1000.0", "");
    EXPECT_THROW(run(), std::exception);
}

TEST_F(Pretension, RejectsNonlinearAnalysis) {
    replace(deck, "TYPE=LINEARSTATIC", "TYPE=NONLINEARSTATIC");
    EXPECT_THROW(run(), std::exception);
}

TEST_F(Pretension, RejectsTopologyAnalysis) {
    replace(deck, "TYPE=LINEARSTATIC", "TYPE=LINEARSTATICTOPO");
    EXPECT_THROW(run(), std::exception);
}

TEST_F(Pretension, ExplicitOutputSelection) {
    replace(deck, "*PRETENSION, SECTION=PT1, ACTION=LOAD", "*NODE FILE\nU\n*PRETENSION, SECTION=PT1, ACTION=LOAD");
    run();
    const auto frames = read_frames(read_file(directory / "model.res"));
    ASSERT_EQ(frames.size(), 3u);
    EXPECT_EQ(frames[0].count("PTGAP"), 0u);
    EXPECT_EQ(frames[0].count("PTFORC"), 0u);
    EXPECT_EQ(frames[0].count("DISPLACEMENT"), 1u);
}

} // namespace

// A larger, obliquely rotated cut checks dense IDs, reference positions and
// every boundary face against its owning element, including the outer faces.
TEST(PretensionPreparation, LargeRotatedCutAndGeneratedResultIds) {
    model::Model model;
    const int n = 16;
    const auto node = [n](int x, int y, int z) { return z*(n+1)*(n+1)+y*(n+1)+x; };
    const Mat3 rotation = Eigen::AngleAxis<Precision>(0.63, Vec3(1,2,3).normalized()).toRotationMatrix();
    for (int z=0; z<4; ++z) for (int y=0; y<=n; ++y) for (int x=0; x<=n; ++x) {
        const Vec3 p = rotation * Vec3(Precision(x)/n, Precision(y)/n, z);
        model.set_node(node(x,y,z), p.x(), p.y(), p.z());
    }
    auto part = model._data->parts.get();
    int e=0, surface=0;
    std::vector<int> face_ids;
    for (int z=0; z<3; ++z) for (int y=0; y<n; ++y) for (int x=0; x<n; ++x) {
        model.set_element<model::C3D8>(e,
            node(x,y,z),node(x+1,y,z),node(x+1,y+1,z),node(x,y+1,z),
            node(x,y,z+1),node(x+1,y,z+1),node(x+1,y+1,z+1),node(x,y+1,z+1));
        part->surface_sets.activate("BOUNDARIES");
        for (int face=1; face<=6; ++face) { model.set_surface(surface++,e,face); face_ids.push_back(face); }
        if (z==1) {
            part->surface_sets.activate("CUT");
            model.set_surface(surface++,e,1);
            face_ids.push_back(1);
        }
        ++e;
    }
    model.compile();
    // A high semantic ID forces batch assignment to scan the mapping rather
    // than assuming that dense and semantic IDs coincide.
    const auto instance = model._data->instances.get(model::Model::DEFAULT_INSTANCE_NAME);
    model._data->node_mapping[0] = {instance, 900001};
    model.add_pretension_interface_section("PT", "CUT", "CUT");
    const auto before = model._data->positions->rows;
    model.prepare_pretension_sections();
    const auto& data = *model._data;
    const auto& section = *data.pretension_sections[0];
    ASSERT_EQ(section.interface_pairs.size(), (n+1)*(n+1));
    EXPECT_EQ(data.positions->rows, before + (n+1)*(n+1) + 1);
    EXPECT_EQ(section.control_node, data.positions->rows-1);
    EXPECT_NEAR(std::abs(section.axis_direction.dot(rotation*Vec3::UnitZ())), 1, 1e-12);
    for (Index i=before; i<data.positions->rows; ++i) {
        EXPECT_EQ(std::get<1>(data.node_mapping[i]), 900002+i-before);
        EXPECT_EQ((data.positions->row_vec3(i)-data.positions_reference->row_vec3(i)).norm(), 0);
    }
    for (const auto& pair : section.interface_pairs) {
        EXPECT_NE(pair.side_a, pair.side_b);
        EXPECT_EQ((data.positions->row_vec3(pair.side_a)-data.positions->row_vec3(pair.side_b)).norm(), 0);
    }
    for (std::size_t i=0; i<data.surfaces.size(); ++i) {
        const auto expected = data.elements[data.surface_element_ids[i]]->surface(face_ids[i]);
        for (Index j=0; j<expected->n_nodes; ++j)
            EXPECT_EQ(data.surfaces[i]->nodes()[j], expected->nodes()[j]);
    }
    model.prepare_pretension_sections();
    EXPECT_EQ(data.positions->rows, before + (n+1)*(n+1) + 1);
}

TEST_F(Pretension, InvalidGeometryDoesNotMutateTopology) {
    replace(deck, "7,   1,  1, 1", "7,   1,  1, 1.2");
    EXPECT_THROW(run(), std::exception);
    const auto& data = *parser.model()._data;
    EXPECT_EQ(data.positions->rows, 16);
    EXPECT_EQ(data.node_mapping.size(), 16);
    EXPECT_FALSE(data.pretension_sections[0]->prepared);
    EXPECT_EQ(data.pretension_sections[0]->control_node, -1);
    // Compiled dense cut nodes are 4..7; no duplicates may be installed.
    for (const auto& element : data.elements)
        for (const ID node : *element) EXPECT_LT(node, 16);
}
