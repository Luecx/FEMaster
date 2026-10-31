#include "../src/io/reader/parser.h"
#include "../src/model/model.h"
#include "../src/bc/structural/pretension.h"
#include "../src/loadcase/linear_static.h"
#include "../src/constraints/transformer/constraint_transformer.h"

#include <gtest/gtest.h>
#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <map>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace fem;

namespace {
std::string read_file(const std::filesystem::path& path) {
    std::ifstream input(path);
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
            row >> value;
            std::vector<double> values;
            while (row >> value) values.push_back(std::strtod(value.c_str(), nullptr));
            frames.back()[field].push_back(std::move(values));
        }
    }
    return frames;
}

void check_section_loadsteps(const std::vector<Frame>& frames,
                            const pretension::PretensionSection& section,
                            const std::vector<double>& expected_forces,
                            const std::vector<double>& expected_gaps = {}) {
    ASSERT_EQ(frames.size(), expected_forces.size());
    ASSERT_TRUE(expected_gaps.empty() || expected_gaps.size() == frames.size());
    const double preload_gap = expected_forces.front() * 3.0 / (210000.0 * 4.0);
    ASSERT_EQ(section.interface_pairs.size(), 4u);
    EXPECT_NEAR(section.axis_direction.norm(), 1, 1e-12);
    for (std::size_t i = 0; i < frames.size(); ++i) {
        SCOPED_TRACE(section.name + " frame " + std::to_string(i));
        const auto& frame = frames[i];
        const double expected_gap = expected_gaps.empty() ? preload_gap : expected_gaps[i];
        double total_force = 0;
        for (const auto& pair : section.interface_pairs) {
            Vec3 relative, force_a, force_b;
            for (Dim c = 0; c < 3; ++c) {
                relative(c) = frame.at("DISPLACEMENT").at(pair.side_b).at(c)
                    - frame.at("DISPLACEMENT").at(pair.side_a).at(c);
                force_a(c) = frame.at("PTFORC").at(pair.side_a).at(c);
                force_b(c) = frame.at("PTFORC").at(pair.side_b).at(c);
            }
            EXPECT_NEAR(section.axis_direction.dot(relative), expected_gap, 1e-8);
            EXPECT_NEAR((relative - expected_gap * section.axis_direction).norm(), 0, 1e-8);
            EXPECT_NEAR((force_a + force_b).norm(), 0, 1e-8);
            for (const ID node : {pair.side_a, pair.side_b})
                EXPECT_NEAR(frame.at("PTGAP").at(node).at(0), expected_gap, 1e-8);
            EXPECT_NEAR(frame.at("TF").at(pair.side_b).at(0), expected_forces[i], 1e-3);
            total_force += force_b.dot(section.axis_direction);
        }
        EXPECT_NEAR(total_force, expected_forces[i], 1e-3);
    }
    EXPECT_NEAR(section.last_solved_gap, expected_gaps.empty() ? preload_gap : expected_gaps.back(), 1e-8);
}

class Pretension : public testing::Test {
protected:
    std::filesystem::path directory;
    std::string deck;
    io::reader::Parser parser;

    void SetUp() override {
        directory = std::filesystem::temp_directory_path() /
            (std::string("femaster_pretension_") + testing::UnitTest::GetInstance()->current_test_info()->name());
        std::filesystem::create_directories(directory);
        deck = read_file(std::filesystem::path(__FILE__).parent_path() / "data/pretension_face_pair_loadsteps.inp");
        replace(deck, "210000.0, 0.3", "210000.0, 0.0");
    }
    void TearDown() override {
        std::filesystem::remove(directory / "model.inp");
        std::filesystem::remove(directory / "model.res");
        std::filesystem::remove(directory / "model.frd");
        std::filesystem::remove(directory / "model.femr");
        std::filesystem::remove(directory);
    }
    void run() {
        std::ofstream(directory / "model.inp") << deck;
        parser.run((directory / "model.inp").string(), (directory / "model").string(), {true, true, false});
    }
    void use_elimination() {
        while (deck.find("TYPE=NULLSPACE") != std::string::npos)
            replace(deck, "TYPE=NULLSPACE", "TYPE=ELIMINATION");
    }
    void use_step_syntax() {
        std::istringstream input(deck);
        std::ostringstream result;
        std::string line;
        while (std::getline(input, line)) {
            if (line.rfind("*LOADCASE, TYPE=LINEARSTATIC, NAME=", 0) == 0)
                result << "*STEP, NAME=" << line.substr(line.find("NAME=") + 5) << "\n*STATIC\n";
            else if (line == "*END") result << "*END STEP\n";
            else result << line << '\n';
        }
        deck = result.str();
    }
    void check_force_hold_force() {
        const auto third = deck.rfind("*PRETENSION, SECTION=PT1, MODE=HOLD");
        ASSERT_NE(third, std::string::npos);
        deck.replace(third, std::string("*PRETENSION, SECTION=PT1, MODE=HOLD").size(),
            "*PRETENSION, SECTION=PT1, MODE=FORCE, VALUE=2000");
        run();
        const double ea = 210000.0 * 4.0;
        check_section_loadsteps(read_frames(read_file(directory / "model.res")),
            *parser.model()._data->pretension_sections.front(),
            {1000, 1000, 2000}, {3000 / ea, 3000 / ea, 5996 / ea});
    }
    void rotate_model() {
        const auto rotation = Eigen::AngleAxis<Precision>(Precision(0.7),
            Vec3(1, 2, 3).normalized()).toRotationMatrix();
        std::istringstream input(deck);
        std::ostringstream rotated;
        rotated << std::setprecision(17);
        std::string line;
        bool nodes = false;
        while (std::getline(input, line)) {
            if (!line.empty() && line.front() == '*') nodes = line.rfind("*NODE,", 0) == 0;
            if ((nodes && !line.empty() && line.front() != '*') || line.rfind("MIDDLE, 0.0, 0.0, -1.0", 0) == 0) {
                std::replace(line.begin(), line.end(), ',', ' ');
                std::istringstream values(line);
                std::string id;
                Vec3 xyz;
                values >> id >> xyz(0) >> xyz(1) >> xyz(2);
                xyz = rotation * xyz;
                rotated << id << ", " << xyz(0) << ", " << xyz(1) << ", " << xyz(2);
                if (!nodes) rotated << ", 0, 0, 0";
                rotated << '\n';
            } else rotated << line << '\n';
        }
        deck = rotated.str();
    }
    void check_single_section() {
        run();
        const auto frames = read_frames(read_file(directory / "model.res"));
        // A load of -4 at z=2 on a length-3 fixed-ended bar increases
        // the compressive cut force at z=1 by 4/3.
        check_section_loadsteps(frames, *parser.model()._data->pretension_sections.front(),
            {1000, 1000, 1000 + 4.0 / 3.0});
    }
    void check_two_sections() {
        run();
        const auto frames = read_frames(read_file(directory / "model.res"));
        const auto& sections = parser.model()._data->pretension_sections;
        ASSERT_EQ(sections.size(), 2u);
        ASSERT_EQ(parser.model()._data->positions->rows, 40);
        check_section_loadsteps(frames, *sections[0], {1000, 1000, 1000 + 4.0 / 3.0});
        check_section_loadsteps(frames, *sections[1], {2000, 2000, 2000});
        EXPECT_NEAR(sections[0]->axis_direction.dot(sections[1]->axis_direction), 0, 1e-12);
    }
    void rejects_cut_constraint() {
        try {
            run();
            FAIL() << "Expected rejection of an additional constraint at the cut";
        } catch (const std::exception& error) {
            EXPECT_NE(std::string(error.what()).find("additional constraint at cut node"),
                std::string::npos) << error.what();
            EXPECT_NE(std::string(error.what()).find("PT1"), std::string::npos);
        }
    }
};

TEST_F(Pretension, ForceHoldExternalLoad) {
    run();
    const auto& data = *parser.model()._data;
    ASSERT_EQ(data.positions->rows, 20);
    ASSERT_EQ(data.pretension_sections.size(), 1u);
    const auto& section = *data.pretension_sections.front();
    ASSERT_EQ(section.interface_pairs.size(), 4u);
    const auto frames = read_frames(read_file(directory / "model.res"));
    ASSERT_EQ(frames.size(), 3u);
    const double expected_gap = 1000.0 * 3.0 / (210000.0 * 4.0);
    for (const auto& frame : frames) {
        const auto& gap = frame.at("PTGAP");
        const auto& displacement = frame.at("DISPLACEMENT");
        for (const auto& pair : section.interface_pairs) {
            EXPECT_NEAR(gap.at(pair.side_a).at(0), expected_gap, 1e-8);
            EXPECT_NEAR(gap.at(pair.side_b).at(0), expected_gap, 1e-8);
            Vec3 relative;
            for (Dim c = 0; c < 3; ++c)
                relative(c) = displacement.at(pair.side_b).at(c) - displacement.at(pair.side_a).at(c);
            EXPECT_NEAR((relative - expected_gap * section.axis_direction).norm(), 0, 1e-8);
        }
    }
    const auto cut_force = [&](const Frame& frame) {
        double total = 0;
        for (const auto& pair : section.interface_pairs)
            for (Dim c = 0; c < 3; ++c)
                total += frame.at("PTFORC").at(pair.side_b).at(c) * section.axis_direction(c);
        return total;
    };
    EXPECT_NEAR(cut_force(frames[0]), 1000.0, 1e-3);
    EXPECT_NEAR(cut_force(frames[1]), 1000.0, 1e-3);
    EXPECT_GT(std::abs(cut_force(frames[2]) - cut_force(frames[1])), 0.1);
}

TEST_F(Pretension, RejectsSharedTopology) {
    replace(deck, "2, 17, 18, 19, 20, 9, 10, 11, 12",
        "2, 5, 6, 7, 8, 9, 10, 11, 12");
    EXPECT_THROW(run(), std::exception);
}
TEST_F(Pretension, RejectsMissingMatch) {
    replace(deck, "17, -1, -1, 1", "17, -1, -1, 1.1");
    EXPECT_THROW(run(), std::exception);
}
TEST_F(Pretension, RejectsAmbiguousMatch) {
    replace(deck, "18,  1, -1, 1", "18, -1, -1, 1");
    EXPECT_THROW(run(), std::exception);
}
TEST_F(Pretension, RejectsHoldWithoutCompletedStep) {
    replace(deck, "*PRETENSION, SECTION=PT1, MODE=FORCE, VALUE=1000.0",
        "*PRETENSION, SECTION=PT1, MODE=HOLD");
    EXPECT_THROW(run(), std::exception);
}
TEST(PretensionPostCheck, RejectsUnbalancedAndNonfiniteSolutions) {
    SystemDofIds dofs(1, 6);
    dofs.setConstant(-1);
    dofs(0, 0) = 0;
    constraint::ConstraintTransformer::Options options;
    options.method = constraint::ConstraintTransformer::Method::NullSpace;
    constraint::ConstraintTransformer transformer({}, dofs, 1, options);
    SparseMatrix stiffness(1, 1);
    stiffness.insert(0, 0) = 1;
    const DynamicVector load = DynamicVector::Ones(1);
    EXPECT_TRUE(transformer.post_check_static(stiffness, load, load));
    EXPECT_FALSE(transformer.post_check_static(stiffness, load, DynamicVector::Zero(1)));
    EXPECT_FALSE(transformer.post_check_static(stiffness, load,
        DynamicVector::Constant(1, std::numeric_limits<Precision>::quiet_NaN())));
}

TEST_F(Pretension, RejectsFailedPostCheckBeforeCapturingGap) {
    replace(deck, "VALUE=1000.0", "VALUE=1e200");
    try {
        run();
        FAIL() << "Expected the numerical post-check to reject the step";
    } catch (const std::exception& error) {
        EXPECT_NE(std::string(error.what()).find("numerical post-check failed"), std::string::npos);
    }
    const auto& section = *parser.model()._data->pretension_sections.front();
    EXPECT_FALSE(section.has_solved_gap);
    EXPECT_EQ(section.last_solved_gap, 0);
    EXPECT_THROW(bc::PretensionCondition(parser.model()._data->pretension_sections.front(),
        bc::PretensionCondition::Mode::Hold), std::exception);
    EXPECT_TRUE(read_frames(read_file(directory / "model.res")).empty());
}
TEST_F(Pretension, FailedStepPreservesPreviouslyAcceptedGap) {
    replace(deck, "*PRETENSION, SECTION=PT1, MODE=HOLD",
        "*PRETENSION, SECTION=PT1, MODE=FORCE, VALUE=1e200");
    try {
        run();
        FAIL() << "Expected the second step to fail its numerical post-check";
    } catch (const std::exception& error) {
        EXPECT_NE(std::string(error.what()).find("numerical post-check failed"), std::string::npos);
    }
    const auto& section = *parser.model()._data->pretension_sections.front();
    const double expected_gap = 1000.0 * 3.0 / (210000.0 * 4.0);
    EXPECT_TRUE(section.has_solved_gap);
    EXPECT_NEAR(section.last_solved_gap, expected_gap, 1e-8);
    const bc::PretensionCondition hold(parser.model()._data->pretension_sections.front(),
        bc::PretensionCondition::Mode::Hold);
    EXPECT_NEAR(hold.value, expected_gap, 1e-8);
    EXPECT_EQ(read_frames(read_file(directory / "model.res")).size(), 1u);
}
TEST_F(Pretension, RejectsActiveSupportOnSideA) {
    replace(deck, "*PRETENSION, SECTION=PT1, MODE=FORCE, VALUE=1000.0",
        "*BOUNDARY\n5, 1, 1, 0\n*PRETENSION, SECTION=PT1, MODE=FORCE, VALUE=1000.0");
    rejects_cut_constraint();
}

TEST_F(Pretension, RejectsActiveSupportOnSideBDuringHold) {
    replace(deck, "*PRETENSION, SECTION=PT1, MODE=HOLD",
        "*BOUNDARY\n17, 1, 1, 0\n*PRETENSION, SECTION=PT1, MODE=HOLD");
    rejects_cut_constraint();
    EXPECT_TRUE(parser.model()._data->pretension_sections.front()->has_solved_gap);
    EXPECT_EQ(read_frames(read_file(directory / "model.res")).size(), 1u);
}

TEST_F(Pretension, RejectsManualEquationAtCut) {
    replace(deck, "*MATERIAL, NAME=STEEL",
        "*EQUATION\n2\n5, 1, 1, 9, 1, -1\n*MATERIAL, NAME=STEEL");
    rejects_cut_constraint();
}

TEST_F(Pretension, RejectsCouplingAtCut) {
    replace(deck, "*MATERIAL, NAME=STEEL",
        "*COUPLING, MASTER=9, SURFACE=CUT_B\n*KINEMATIC\n1, 3\n*MATERIAL, NAME=STEEL");
    rejects_cut_constraint();
}

TEST_F(Pretension, AllowsInactiveSupportDefinitionAtCut) {
    replace(deck, "*MATERIAL, NAME=STEEL",
        "*NSET, NAME=CUT_NODES\n5, 17\n*MATERIAL, NAME=STEEL");
    replace(deck, "*CLOAD, LOAD_COLLECTOR=EXTERNAL",
        "*SUPPORT, SUPPORT_COLLECTOR=UNUSED\nCUT_NODES, 0, 0, 0, NAN, NAN, NAN\n"
        "*CLOAD, LOAD_COLLECTOR=EXTERNAL");
    EXPECT_NO_THROW(run());
}

TEST_F(Pretension, AllowsExternalLoadAtCutAndRecoversForce) {
    replace(deck, "*CLOAD, LOAD_COLLECTOR=EXTERNAL",
        "*NSET, NAME=LOADED_CUT\n17, 18, 19, 20\n*CLOAD, LOAD_COLLECTOR=CUT_LOAD\n"
        "LOADED_CUT, 0, 0, 10, 0, 0, 0\n*CLOAD, LOAD_COLLECTOR=EXTERNAL");
    replace(deck, "*PRETENSION, SECTION=PT1, MODE=FORCE, VALUE=1000.0",
        "*LOADS\nCUT_LOAD\n*PRETENSION, SECTION=PT1, MODE=FORCE, VALUE=1000.0");
    run();
    const auto frames = read_frames(read_file(directory / "model.res"));
    ASSERT_EQ(frames.size(), 3u);
    const auto& section = *parser.model()._data->pretension_sections.front();
    Precision force = 0;
    for (const auto& pair : section.interface_pairs)
        for (Dim c = 0; c < 3; ++c)
            force += frames[0].at("PTFORC").at(pair.side_b).at(c) * section.axis_direction(c);
    EXPECT_NEAR(force, 1000, 1e-3);
}
TEST_F(Pretension, AxialForceHoldWithElimination) {
    use_elimination();
    check_single_section();
}

TEST_F(Pretension, ObliqueForceHoldWithNullspace) {
    rotate_model();
    check_single_section();
}

TEST_F(Pretension, ObliqueForceHoldWithElimination) {
    rotate_model();
    use_elimination();
    check_single_section();
}

TEST_F(Pretension, TwoIndependentSectionsWithNullspace) {
    deck = read_file(std::filesystem::path(__FILE__).parent_path() / "data/pretension_two_sections.inp");
    check_two_sections();
}

TEST_F(Pretension, TwoIndependentSectionsWithElimination) {
    deck = read_file(std::filesystem::path(__FILE__).parent_path() / "data/pretension_two_sections.inp");
    use_elimination();
    check_two_sections();
}

TEST_F(Pretension, HeldGapWithChangedExternalLoad) {
    replace(deck, "MIDDLE, 0.0, 0.0, -1.0", "MIDDLE, 0.0, 0.0, -2.0");
    run();
    const auto frames = read_frames(read_file(directory / "model.res"));
    check_section_loadsteps(frames, *parser.model()._data->pretension_sections.front(),
        {1000, 1000, 1000 + 8.0 / 3.0});
}
TEST_F(Pretension, ReplacementPreservesFrozenInitialPrescription) {
    run();
    using Prescription = bc::PretensionCondition;
    const auto section = parser.model()._data->pretension_sections.front();
    const auto accepted = *parser.conditions().get(bc::PRETENSION).begin();
    const auto& held = static_cast<const Prescription&>(*accepted);
    ASSERT_EQ(held.mode, Prescription::Mode::Hold);
    const auto gap = held.value;
    parser.begin_loadcase(std::make_unique<loadcase::LinearStatic>());
    const auto force = std::make_shared<Prescription>(section, Prescription::Mode::Force, 2500);
    parser.modify_conditions(bc::PRETENSION, section->name, {force});
    const auto* step = parser.active_loadcase();
    EXPECT_TRUE(step->initial_conditions.contains(bc::PRETENSION, accepted));
    EXPECT_FALSE(step->initial_conditions.contains(bc::PRETENSION, force));
    EXPECT_FALSE(step->target_conditions.contains(bc::PRETENSION, accepted));
    EXPECT_TRUE(step->target_conditions.contains(bc::PRETENSION, force));
    EXPECT_EQ(held.value, gap);
    const auto hold = std::make_shared<Prescription>(section, Prescription::Mode::Hold);
    parser.modify_conditions(bc::PRETENSION, section->name, {hold});
    ASSERT_EQ(step->target_conditions.get(bc::PRETENSION).size(), 1u);
    EXPECT_TRUE(step->initial_conditions.contains(bc::PRETENSION, accepted));
    EXPECT_TRUE(step->target_conditions.contains(bc::PRETENSION, hold));
    EXPECT_EQ(hold->value, gap);
}

TEST_F(Pretension, RejectedTargetDoesNotBecomeInheritedPrescription) {
    run();
    using Prescription = bc::PretensionCondition;
    const auto section = parser.model()._data->pretension_sections.front();
    const Precision gap = section->last_solved_gap;
    const auto accepted = *parser.conditions().get(bc::PRETENSION).begin();
    auto trial = std::make_unique<loadcase::LinearStatic>();
    trial->constraint_method = constraint::ConstraintTransformer::Method::NullSpace;
    parser.begin_loadcase(std::move(trial));
    auto* step = parser.active_loadcase();
    for (const auto& support : step->initial_conditions.get(bc::SUPPORT))
        step->target_conditions.add(bc::SUPPORT, support);
    const auto rejected = std::make_shared<Prescription>(section, Prescription::Mode::Force, 1e200);
    parser.modify_conditions(bc::PRETENSION, section->name, {rejected});
    try {
        parser.end_loadcase();
        FAIL() << "Expected a rejected numerical solution";
    } catch (const std::exception& error) {
        EXPECT_NE(std::string(error.what()).find("numerical post-check failed"), std::string::npos);
    }
    EXPECT_EQ(section->last_solved_gap, gap);
    EXPECT_TRUE(parser.conditions().contains(bc::PRETENSION, accepted));
    EXPECT_FALSE(parser.conditions().contains(bc::PRETENSION, rejected));
    parser.begin_loadcase(std::make_unique<loadcase::LinearStatic>());
    EXPECT_TRUE(parser.active_loadcase()->initial_conditions.contains(bc::PRETENSION, accepted));
    EXPECT_TRUE(parser.active_loadcase()->target_conditions.contains(bc::PRETENSION, accepted));
    EXPECT_FALSE(parser.active_loadcase()->target_conditions.contains(bc::PRETENSION, rejected));
}

TEST_F(Pretension, InheritsForceWithoutRepeatedPrescription) {
    while (deck.find("*PRETENSION, SECTION=PT1, MODE=HOLD\n") != std::string::npos)
        replace(deck, "*PRETENSION, SECTION=PT1, MODE=HOLD\n", "");
    run();
    const double ea = 210000.0 * 4.0;
    check_section_loadsteps(read_frames(read_file(directory / "model.res")),
        *parser.model()._data->pretension_sections.front(),
        {1000, 1000, 1000}, {3000 / ea, 3000 / ea, 2996 / ea});
}

TEST_F(Pretension, InheritsHoldWithoutRepeatedPrescription) {
    const auto third = deck.rfind("*PRETENSION, SECTION=PT1, MODE=HOLD\n");
    ASSERT_NE(third, std::string::npos);
    deck.erase(third, std::string("*PRETENSION, SECTION=PT1, MODE=HOLD\n").size());
    check_single_section();
}

TEST_F(Pretension, ForceHoldForceWithLoadcaseSyntax) {
    check_force_hold_force();
}

TEST_F(Pretension, ForceReplacementDoesNotAccumulate) {
    replace(deck, "*PRETENSION, SECTION=PT1, MODE=HOLD",
        "*PRETENSION, SECTION=PT1, MODE=FORCE, VALUE=2000");
    run();
    const double ea = 210000.0 * 4.0;
    check_section_loadsteps(read_frames(read_file(directory / "model.res")),
        *parser.model()._data->pretension_sections.front(),
        {1000, 2000, 2000 + 4.0 / 3.0}, {3000 / ea, 6000 / ea, 6000 / ea});
}

TEST_F(Pretension, ForceHoldForceWithStepSyntax) {
    use_step_syntax();
    check_force_hold_force();
}
} // namespace
