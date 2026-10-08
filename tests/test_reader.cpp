/**
 * @file test_reader.cpp
 * @brief Verifies parser registration and selected material keyword mappings.
 *
 * The tests exercise rigid-body-motion command registration and parsing,
 * orthotropic engineering-constant mapping, and shell shear-component ordering
 * through the public reader and material interfaces. Temporary input and result
 * files are removed before and after each parser scenario.
 *
 * @see io::reader::Parser
 * @see material::OrthotropicElasticity
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#include "../src/io/reader/parser.h"
#include "../src/io/dsl/deck_parser.h"
#include "../src/bc/structural/load_c.h"
#include "../src/bc/amplitude.h"
#include "../src/bc/structural/load_p.h"
#include "../src/bc/structural/support.h"
#include "../src/bc/structural/load_inertial.h"
#include "../src/bc/thermal/temperature.h"
#include "../src/loadcase/linear_static.h"
#include "../src/loadcase/nonlinear_static.h"
#include "../src/material/orthotropic_elasticity.h"
#include "../src/material/strain/shell_material_strain_green_lagrange.h"
#include "../src/material/stress/shell_material_stress_pk2.h"
#include "../src/model/model.h"
#include "../src/section/section_solid.h"
#include "../src/model/solid/c3d5.h"
#include "../src/section/section_truss.h"
#include "../src/section/section_shell_integrated.h"

#include <filesystem>
#include <fstream>
#include <cmath>
#include <vector>
#include <memory>
#include <stdexcept>
#include <utility>
#include <string>

#include <gtest/gtest.h>

using namespace fem;

TEST(Reader_Parser, RegistersRbmCommand) {
    io::reader::Parser parser;
    EXPECT_NE(parser.registry().find("RBM"), nullptr);
}


TEST(Reader_Parser, NativeStepNestsExactlyOneProcedureAndItsConditions) {
    const std::string path = "tests/TMP_NATIVE_STEP_SCOPE.INP";
    {
        std::ofstream out(path);
        ASSERT_TRUE(out.is_open());
        out << "*STEP, NAME=TEST, NLGEOM=YES\n";
        out << "*STATIC\n";
        out << "*CLOAD\n";
        out << "NODES, 1, 15.0\n";
        out << "*BOUNDARY\n";
        out << "NODES, 1, 3\n";
        out << "*OUTPUT, FIELD\n";
        out << "*NODE OUTPUT\n";
        out << "U\n";
        out << "*END STEP\n";
    }

    io::reader::Parser parser;
    io::dsl::File file(path);
    io::dsl::DeckParser deck_parser(parser.registry());
    const auto deck = deck_parser.parse(file);
    const auto steps = deck.root().children("STEP");
    ASSERT_EQ(steps.size(), 1u);

    const auto procedures = steps.front()->children();
    ASSERT_EQ(procedures.size(), 1u);
    EXPECT_EQ(procedures.front()->command().name_, "STATIC");

    const auto commands = procedures.front()->children();
    ASSERT_EQ(commands.size(), 5u);
    EXPECT_EQ(commands[0]->command().name_, "CLOAD");
    EXPECT_EQ(commands[1]->command().name_, "BOUNDARY");
    EXPECT_EQ(commands[2]->command().name_, "OUTPUT");
    EXPECT_EQ(commands[3]->command().name_, "NODEOUTPUT");
    EXPECT_EQ(commands[4]->command().name_, "ENDSTEP");

    std::filesystem::remove(path);
}

TEST(Reader_Parser, RejectsStepAmplitudeModes) {
    const std::string path = "tests/TMP_STEP_AMPLITUDE.INP";
    {
        std::ofstream out(path);
        ASSERT_TRUE(out.is_open());
        out << "*STEP, AMPLITUDE=RAMP\n";
        out << "*STATIC\n";
        out << "*END STEP\n";
        out << "*STEP, AMPLITUDE=STEP\n";
        out << "*STATIC\n";
        out << "*END STEP\n";
    }

    io::reader::Parser parser;
    io::dsl::File file(path);
    io::dsl::DeckParser deck_parser(parser.registry());
    const auto deck = deck_parser.parse(file);

    const auto steps = deck.root().children("STEP");
    ASSERT_EQ(steps.size(), 2u);
    for (const auto* step : steps) {
        EXPECT_THROW(step->enter(), std::exception);
    }

    std::filesystem::remove(path);
}

TEST(Reader_Parser, LinearStaticSamplesNamedAmplitudeAtStepEnd) {
    const std::string path = "tests/TMP_LINEAR_STATIC_AMPLITUDE.INP";
    {
        std::ofstream out(path);
        ASSERT_TRUE(out.is_open());
        out << "*AMPLITUDE, NAME=CURVE\n";
        out << "0.0, 0.0, 2.5, 1.0\n";
        out << "*STEP\n*STATIC\n0.25, 2.5\n*END STEP\n";
    }

    io::reader::Parser parser;
    {
        io::dsl::File file(path);
        io::dsl::DeckParser deck_parser(parser.registry());
        const auto deck = deck_parser.parse(file);

        parser.model().set_node(1, Precision(0), Precision(0), Precision(0));
        parser.model().compile();
        deck.root().execute_children("AMPLITUDE");

        const auto steps = deck.root().children("STEP");
        ASSERT_EQ(steps.size(), 1u);
        const auto procedures = steps.front()->children("STATIC");
        ASSERT_EQ(procedures.size(), 1u);

        steps.front()->enter();
        procedures.front()->enter();

        auto* active = parser.active_loadcase();
        ASSERT_NE(active, nullptr);
        auto* linear = active->as<loadcase::LinearStatic>();
        ASSERT_NE(linear, nullptr);
        EXPECT_DOUBLE_EQ(linear->step_period, 2.5);

        auto load = std::make_shared<bc::CLoad>();
        load->region_     = parser.model().resolve_node_region("1");
        load->values_[2]  = Precision(-1000);
        load->amplitude_  = parser.model()._data->amplitudes.get("CURVE");
        parser.add_condition(bc::CLOAD, load, "NODE:1");

        const auto at_start = parser.model().build_load_matrix(Precision(0));
        const auto at_end   = parser.model().build_load_matrix(linear->step_period);
        const ID   node_id  = parser.model().compiled_node_id(1);

        EXPECT_DOUBLE_EQ(at_start(node_id, 2), Precision(0));
        EXPECT_DOUBLE_EQ(at_end(node_id, 2), Precision(-1000));
    }

    std::filesystem::remove(path);
}

TEST(Reader_Parser, NativeReaderAcceptsAbaqusReducedShellAliases) {
    const std::string path = "tests/TMP_SHELL_ALIASES.INP";
    {
        std::ofstream out(path);
        ASSERT_TRUE(out.is_open());
        out << "*ELEMENT, TYPE=S3R\n1, 1, 2, 3\n";
        out << "*ELEMENT, TYPE=S4R\n2, 1, 2, 3, 4\n";
        out << "*ELEMENT, TYPE=S6R\n3, 1, 2, 3, 4, 5, 6\n";
        out << "*ELEMENT, TYPE=S8R\n4, 1, 2, 3, 4, 5, 6, 7, 8\n";
    }

    io::reader::Parser parser;
    io::dsl::File file(path);
    io::dsl::DeckParser deck_parser(parser.registry());
    const auto deck = deck_parser.parse(file);
    EXPECT_EQ(deck.root().children("ELEMENT").size(), 4u);

    std::filesystem::remove(path);
}

TEST(Reader_Parser, C3D5CreatesDedicatedPyramidElement) {
    const std::string path = "tests/TMP_C3D5_ELEMENT.INP";
    {
        std::ofstream out(path);
        ASSERT_TRUE(out.is_open());
        out << "*NODE\n"
            << "1, 0, 0, 0\n"
            << "2, 1, 0, 0\n"
            << "3, 1, 1, 0\n"
            << "4, 0, 1, 0\n"
            << "5, 0.5, 0.5, 1\n"
            << "*ELEMENT, TYPE=C3D5, ELSET=PYRAMID\n"
            << "10, 1, 2, 3, 4, 5\n";
    }

    io::reader::Parser parser;
    io::dsl::File file(path);
    io::dsl::DeckParser deck_parser(parser.registry());
    const auto deck = deck_parser.parse(file);

    deck.root().execute_children("NODE");
    deck.root().execute_children("ELEMENT");

    const auto part = parser.model()._data->parts.get();
    ASSERT_NE(part, nullptr);
    const auto element = part->elements.at(10);
    ASSERT_NE(element->as<model::C3D5>(), nullptr);
    EXPECT_EQ(element->type_name(), "C3D5");

    std::filesystem::remove(path);
}

TEST(Reader_Parser, OrientationAcceptsTypeOrSystemButNotBoth) {
    const std::string input_path  = "tests/TMP_ORIENTATION_SYSTEM.INP";
    const std::string output_path = "tests/TMP_ORIENTATION_SYSTEM.RES";

    {
        std::ofstream out(input_path);
        ASSERT_TRUE(out.is_open());
        out << "*ORIENTATION, NAME=NATIVE, TYPE=RECTANGULAR\n1, 0, 0\n";
        out << "*ORIENTATION, NAME=ABAQUS, SYSTEM=RECTANGULAR, DEFINITION=COORDINATES\n";
        out << "1, 0, 0, 0, 1, 0\n";
        out << "*ORIENTATION, NAME=DEFAULT\n";
        out << "0, 1, 0, -1, 0, 0\n";
        out << "*ORIENTATION, NAME=COORDINATES, DEFINITION=COORDINATES\n";
        out << "2, 1, 0, 1, 2, 0, 1, 1, 0\n";
    }

    io::reader::Parser parser;
    ASSERT_NO_THROW(parser.run(input_path, output_path));
    EXPECT_TRUE(parser.model()._data->coordinate_systems.has("NATIVE"));
    EXPECT_TRUE(parser.model()._data->coordinate_systems.has("ABAQUS"));
    ASSERT_TRUE(parser.model()._data->coordinate_systems.has("DEFAULT"));
    ASSERT_TRUE(parser.model()._data->coordinate_systems.has("COORDINATES"));

    const auto default_axes = parser.model()._data->coordinate_systems.get("DEFAULT")->get_axes(Vec3::Zero());
    EXPECT_NEAR(default_axes(0, 0), Precision(0), 1e-12);
    EXPECT_NEAR(default_axes(1, 0), Precision(1), 1e-12);

    const auto coordinate_axes = parser.model()._data->coordinate_systems.get("COORDINATES")->get_axes(Vec3::Zero());
    EXPECT_NEAR(coordinate_axes(0, 0), Precision(1), 1e-12);
    EXPECT_NEAR(coordinate_axes(1, 0), Precision(0), 1e-12);

    {
        std::ofstream out(input_path);
        ASSERT_TRUE(out.is_open());
        out << "*ORIENTATION, NAME=INVALID, TYPE=RECTANGULAR, SYSTEM=RECTANGULAR\n";
        out << "1, 0, 0, 0, 1, 0\n";
    }
    EXPECT_THROW(parser.run(input_path, output_path), std::exception);

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);
}

TEST(Reader_Parser, RiksRejectsNamedLoadAmplitudes) {
    io::reader::Parser parser;

    auto nonlinear = std::make_unique<loadcase::NonlinearStatic>();
    nonlinear->control = loadcase::NonlinearControl::ArcLength;
    parser.begin_loadcase(std::move(nonlinear));

    auto load = std::make_shared<bc::CLoad>();
    load->amplitude_ = std::make_shared<bc::Amplitude>("HISTORY");
    parser.add_condition(bc::CLOAD, std::move(load), "NODE:1");

    ASSERT_NE(parser.active_loadcase(), nullptr);
    EXPECT_THROW(parser.active_loadcase()->run(), std::exception);
}

TEST(Reader_Parser, BoundaryCreatesDirectStructuralSupportsAndRejectsThermalDofs) {
    const std::string input_path  = "tests/TMP_BOUNDARY_SUPPORT.INP";
    const std::string output_path = "tests/TMP_BOUNDARY_SUPPORT.RES";

    {
        std::ofstream out(input_path);
        ASSERT_TRUE(out.is_open());
        out << "*NODE\n1, 0, 0, 0\n";
        out << "*BOUNDARY\n1, 1, 3\n";
    }

    io::reader::Parser parser;
    ASSERT_NO_THROW(parser.run(input_path, output_path));
    EXPECT_EQ(parser.conditions().get(bc::SUPPORT).size(), 3u);

    {
        std::ofstream out(input_path);
        ASSERT_TRUE(out.is_open());
        out << "*NODE\n1, 0, 0, 0\n";
        out << "*BOUNDARY\n1, 11, 11\n";
    }
    EXPECT_THROW(parser.run(input_path, output_path), std::exception);

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);
}

TEST(Reader_Parser, MixedSolidSectionCreatesSolidAndTrussSubsets) {
    const std::string input_path  = "tests/TMP_MIXED_SOLID_SECTION.INP";
    const std::string output_path = "tests/TMP_MIXED_SOLID_SECTION.RES";

    {
        std::ofstream out(input_path);
        ASSERT_TRUE(out.is_open());
        out << "*MATERIAL, NAME=MAT\n*ELASTIC\n210000., 0.3\n";
        out << "*NODE\n1, 0, 0, 0\n2, 1, 0, 0\n3, 0, 1, 0\n4, 0, 0, 1\n";
        out << "*ELEMENT, TYPE=C3D4\n1, 1, 2, 3, 4\n";
        out << "*ELEMENT, TYPE=T3\n2, 1, 2\n";
        out << "*ELSET, ELSET=MIXED\n1, 2\n";
        out << "*SOLID SECTION, ELSET=MIXED, MATERIAL=MAT\n12.5\n";
    }

    io::reader::Parser parser;
    ASSERT_NO_THROW(parser.run(input_path, output_path));

    const auto part = parser.model()._data->parts.get();
    ASSERT_NE(part, nullptr);
    ASSERT_EQ(part->sections.size(), 2u);
    ASSERT_NE(part->sections[0]->as<SolidSection>(), nullptr);
    ASSERT_NE(part->sections[1]->as<TrussSection>(), nullptr);
    EXPECT_EQ(part->sections[0]->region_->size(), 1u);
    EXPECT_EQ(part->sections[1]->region_->size(), 1u);
    EXPECT_EQ(part->sections[0]->region_->at(0), 1);
    EXPECT_EQ(part->sections[1]->region_->at(0), 2);
    EXPECT_NE(part->sections[0]->region_, part->sections[1]->region_);
    EXPECT_DOUBLE_EQ(part->sections[1]->as<TrussSection>()->area_, 12.5);
    EXPECT_EQ(part->elem_sets.get("MIXED")->size(), 2u);

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);
}

TEST(Reader_Parser, SolidSectionRequiresAreaOnlyForTrusses) {
    const std::string input_path  = "tests/TMP_SOLID_SECTION_AREA.INP";
    const std::string output_path = "tests/TMP_SOLID_SECTION_AREA.RES";
    io::reader::Parser parser;

    {
        std::ofstream out(input_path);
        ASSERT_TRUE(out.is_open());
        out << "*MATERIAL, NAME=MAT\n*ELASTIC\n210000., 0.3\n";
        out << "*NODE\n1, 0, 0, 0\n2, 1, 0, 0\n3, 0, 1, 0\n4, 0, 0, 1\n";
        out << "*ELEMENT, TYPE=C3D4\n1, 1, 2, 3, 4\n";
        out << "*ELSET, ELSET=SOLIDS\n1\n";
        out << "*SOLID SECTION, ELSET=SOLIDS, MATERIAL=MAT\n";
    }
    ASSERT_NO_THROW(parser.run(input_path, output_path));
    ASSERT_EQ(parser.model()._data->parts.get()->sections.size(), 1u);
    EXPECT_NE(parser.model()._data->parts.get()->sections[0]->as<SolidSection>(), nullptr);

    {
        std::ofstream out(input_path);
        ASSERT_TRUE(out.is_open());
        out << "*MATERIAL, NAME=MAT\n*ELASTIC\n210000., 0.3\n";
        out << "*NODE\n1, 0, 0, 0\n2, 1, 0, 0\n";
        out << "*ELEMENT, TYPE=T3\n1, 1, 2\n";
        out << "*ELSET, ELSET=TRUSSES\n1\n";
        out << "*SOLID SECTION, ELSET=TRUSSES, MATERIAL=MAT\n";
    }
    EXPECT_THROW(parser.run(input_path, output_path), std::exception);

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);
}

TEST(Reader_Parser, ShellSectionRequiresFiveIntegrationPoints) {
    const std::string input_path  = "tests/TMP_SHELL_SECTION_IP.INP";
    const std::string output_path = "tests/TMP_SHELL_SECTION_IP.RES";
    io::reader::Parser parser;

    const auto write_input = [&](int integration_points) {
        std::ofstream out(input_path);
        out << "*MATERIAL, NAME=MAT\n*ELASTIC\n210000., 0.3\n";
        out << "*NODE\n1, 0, 0, 0\n2, 1, 0, 0\n3, 0, 1, 0\n";
        out << "*ELEMENT, TYPE=S3\n1, 1, 2, 3\n";
        out << "*ELSET, ELSET=SHELLS\n1\n";
        out << "*SHELL SECTION, ELSET=SHELLS, MATERIAL=MAT\n2.0, "
            << integration_points << "\n";
    };

    write_input(5);
    ASSERT_NO_THROW(parser.run(input_path, output_path));
    ASSERT_EQ(parser.model()._data->parts.get()->sections.size(), 1u);
    EXPECT_NE(parser.model()._data->parts.get()->sections[0]->as<IntegratedShellSection>(), nullptr);

    write_input(3);
    EXPECT_THROW(parser.run(input_path, output_path), std::exception);

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);
}

TEST(Reader_Parser, RegistersInertialLoadFromNativeDeck) {
    const std::string input_path = "tests/TMP_INERTIALOAD.INP";
    const std::string output_path = "tests/TMP_INERTIALOAD.RES";

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);

    {
        std::ofstream os(input_path);
        ASSERT_TRUE(os.is_open());
        os << "*NODE\n";
        os << "1, 0.0, 0.0, 0.0\n";
        os << "*INERTIALOAD, LOAD_COLLECTOR=GRAVITY\n";
        os << "EALL, 0., 0., 0., 0., 0., 9.81, 0., 0., 0., 0., 0., 0.\n";
    }

    io::reader::Parser parser;
    ASSERT_NO_THROW(parser.run(input_path, output_path));

    auto collector = parser.model()._data->load_cols.get("GRAVITY");
    ASSERT_NE(collector, nullptr);
    ASSERT_EQ(collector->size(), 1u);

    auto load = std::dynamic_pointer_cast<bc::InertialLoad>(collector->first());
    ASSERT_NE(load, nullptr);
    EXPECT_EQ(load->region_, parser.model()._data->elem_sets.get("EALL"));
    EXPECT_DOUBLE_EQ(load->center_acc_(2), 9.81);

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);
}

TEST(Reader_Conditions, RejectsUnsupportedFollowerLoads) {
    const std::string input_path  = "tests/TMP_UNSUPPORTED_FOLLOWER.INP";
    const std::string output_path = "tests/TMP_UNSUPPORTED_FOLLOWER.RES";

    const auto rejects = [&](const std::string& input, const std::string& reason) {
        std::filesystem::remove(input_path);
        std::filesystem::remove(output_path);
        {
            std::ofstream os(input_path);
            ASSERT_TRUE(os.is_open());
            os << "*NODE\n1, 0., 0., 0.\n" << input;
        }

        io::reader::Parser parser;
        try {
            parser.run(input_path, output_path);
            FAIL() << "Expected the reader to reject an unsupported follower load";
        } catch (const std::runtime_error& error) {
            EXPECT_NE(std::string(error.what()).find(reason), std::string::npos)
                << "Unexpected reader error: " << error.what();
        }

        std::filesystem::remove(input_path);
        std::filesystem::remove(output_path);
    };

    // Explicit FOLLOWER=YES is invalid even outside a nonlinear step.
    rejects("*DLOAD, NAME=TEST, FOLLOWER=YES\nMISSING, TRVEC1, 1., 1., 0., 0.\n",
            "FOLLOWER");
    rejects("*DSLOAD, NAME=TEST, FOLLOWER=YES\nMISSING, TRVEC, 1., 1., 0., 0.\n",
            "FOLLOWER");

    // TRVEC is rejected for all procedures, including collector definitions.
    rejects("*DLOAD, NAME=TEST\nMISSING, TRVEC1, 1., 1., 0., 0.\n",
            "TRVEC is not supported");
    rejects("*DSLOAD, NAME=TEST\nMISSING, TRVEC, 1., 1., 0., 0.\n",
            "TRVEC is not supported");

    // Native pressure remains unsupported inside nonlinear analysis.
    rejects("*LOADCASE, TYPE=NONLINEARSTATIC\n*PLOAD\nMISSING, 1.\n*END\n",
            "follower pressure is not supported");
}

TEST(Reader_Conditions, RejectsCollectorPressureInNonlinearStep) {
    io::reader::Parser parser;
    parser.begin_loadcase(std::make_unique<loadcase::NonlinearStatic>());

    auto pressure = std::make_shared<bc::PLoad>();
    EXPECT_THROW(parser.select_collector_condition(bc::CLOAD, pressure), std::runtime_error);

    // Fixed-direction tractions do not require a follower load stiffness and
    // must remain selectable (the normal collector lifecycle is unchanged).
    auto traction = std::make_shared<bc::CLoad>();
    EXPECT_NO_THROW(parser.select_collector_condition(bc::CLOAD, traction));
}

TEST(Reader_Conditions, NativeSupportSplitsNodalTransformsAndReplacesAllFragments) {
    const std::string input_path  = "tests/TMP_SUPPORT_TRANSFORM.INP";
    const std::string output_path = "tests/TMP_SUPPORT_TRANSFORM.RES";

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);

    {
        std::ofstream os(input_path);
        ASSERT_TRUE(os.is_open());
        os << "*NODE\n";
        os << "1, 0., 0., 0.\n";
        os << "2, 1000., 0., 0.\n";
        os << "*NSET, NAME=BOTH\n1, 2\n";
        os << "*NSET, NAME=ROTATED\n2\n";
        os << "*TRANSFORM, NSET=ROTATED, TYPE=R\n";
        os << "0., 1., 0., -1., 0., 0.\n";
        os << "*SUPPORT, SUPPORT_COLLECTOR=BC\n";
        os << "BOTH, 0., 0.\n";
    }

    io::reader::Parser parser;
    ASSERT_NO_THROW(parser.run(input_path, output_path));

    const auto collector = parser.model()._data->supp_cols.get("BC");
    ASSERT_NE(collector, nullptr);
    ASSERT_EQ(collector->size(), 4u);

    // The two prescribed DOFs each generate one global and one transformed
    // fragment. No support may constrain both nodes in a single basis.
    int global_count = 0;
    int local_count  = 0;
    const auto region = parser.model()._data->node_sets.get("BOTH");
    auto& model_data = *parser.model()._data;
    auto rhs = *model_data.positions;

    for (const auto& condition : *collector) {
        const auto support = std::dynamic_pointer_cast<bc::Support>(condition);
        ASSERT_NE(support, nullptr);
        ASSERT_NE(support->node_region(), nullptr);
        ASSERT_EQ(support->node_region()->size(), 1u);

        // Check the actual constraint equation, not just the stored basis
        // name. The local x axis is global +Y; the local y axis is global -X.
        constraint::Equations equations;
        SystemDofIds dof_ids;
        TripletList triplets;
        support->apply(model_data, rhs, equations, dof_ids, triplets, Precision(0));
        ASSERT_EQ(equations.size(), 1u);

        const bool local_x = std::isfinite(support->values()[0]);
        if (support->str().find("orientation=__TRANSFORM_ROTATED") != std::string::npos) {
            ++local_count;
            ASSERT_EQ(equations.front().entries.size(), 3u);
            const auto& entries = equations.front().entries;
            EXPECT_NEAR(entries[0].coeff, local_x ? 0.0 : -1.0, 1e-12);
            EXPECT_NEAR(entries[1].coeff, local_x ? 1.0 : 0.0, 1e-12);
            EXPECT_NEAR(entries[2].coeff, 0.0, 1e-12);
        } else {
            ++global_count;
            EXPECT_EQ(support->str().find("orientation="), std::string::npos);
            ASSERT_EQ(equations.front().entries.size(), 1u);
            EXPECT_EQ(equations.front().entries[0].dof, local_x ? 0 : 1);
            EXPECT_DOUBLE_EQ(equations.front().entries[0].coeff, 1.0);
        }

        // Exercise the same original-source index used by direct SUPPORT
        // conditions. Both transformed fragments must be replaced together.
        parser.add_condition(bc::SUPPORT, condition, "NSET:BOTH");
    }

    EXPECT_EQ(global_count, 2);
    EXPECT_EQ(local_count, 2);

    Vec6 replacement_values = Vec6::Constant(NAN);
    replacement_values[0] = Precision(1);
    parser.modify_conditions(bc::SUPPORT, "NSET:BOTH", {
        std::make_shared<bc::Support>(region, replacement_values)
    });

    // Both old x fragments are replaced by one new x definition; y remains
    // independently active on its original two transformed fragments.
    ASSERT_EQ(parser.conditions().get(bc::SUPPORT).size(), 3u);
    int x_count = 0;
    int y_count = 0;
    for (const auto& condition : parser.conditions().get(bc::SUPPORT)) {
        const auto support = std::dynamic_pointer_cast<bc::Support>(condition);
        ASSERT_NE(support, nullptr);
        x_count += std::isfinite(support->values()[0]) ? 1 : 0;
        y_count += std::isfinite(support->values()[1]) ? 1 : 0;
    }
    EXPECT_EQ(x_count, 1);
    EXPECT_EQ(y_count, 2);

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);
}

TEST(Reader_Conditions, NewPreservesLoadsDefinedInCurrentStep) {
    io::reader::Parser parser;

    auto inherited = std::make_shared<bc::CLoad>();
    inherited->values_[0] = Precision(100);
    parser.add_condition(bc::CLOAD, inherited, "NODE:1");

    parser.begin_loadcase(std::make_unique<loadcase::LinearStatic>());

    auto current = std::make_shared<bc::CLoad>();
    current->values_[0] = Precision(200);
    parser.add_condition(bc::CLOAD, current, "NODE:2");

    parser.clear_conditions(bc::CLOAD);

    // NEW retires the definition inherited at step entry, not the new load.
    EXPECT_TRUE(parser.conditions().contains(bc::CLOAD, inherited));
    EXPECT_TRUE(parser.conditions().contains(bc::CLOAD, current));
    EXPECT_DOUBLE_EQ(inherited->values_[0], Precision(0));
    EXPECT_DOUBLE_EQ(current->values_[0], Precision(200));

    // The current-step source remains indexed and survives another MOD.
    auto replacement = std::make_shared<bc::CLoad>();
    replacement->values_[0] = Precision(300);
    parser.modify_conditions(bc::CLOAD, "NODE:2", {replacement});
    EXPECT_DOUBLE_EQ(current->values_[0], Precision(0));
    EXPECT_DOUBLE_EQ(replacement->values_[0], Precision(300));
}

TEST(Reader_Conditions, ModMutatesCompatibleLoadAndPreservesStepStart) {
    io::reader::Parser parser;
    auto region = std::make_shared<model::NodeRegion>("SOURCE");
    region->add(0);

    auto original = std::make_shared<bc::CLoad>();
    original->region_ = region;
    original->values_[0] = Precision(100);
    parser.add_condition(bc::CLOAD, original, "NSET:SOURCE");

    parser.begin_loadcase(std::make_unique<loadcase::LinearStatic>());
    // Represent an already applied load inherited from a completed step.
    original->values_start_[0] = Precision(100);
    EXPECT_DOUBLE_EQ(original->values_start_[0], Precision(100));

    for (const Precision magnitude : {Precision(200), Precision(300)}) {
        auto incoming = std::make_shared<bc::CLoad>();
        incoming->region_ = region;
        incoming->values_[0] = magnitude;
        parser.modify_conditions(bc::CLOAD, "NSET:SOURCE", {incoming});

        // No duplicate object: only its target changes.
        EXPECT_EQ(parser.conditions().get(bc::CLOAD).size(), 1u);
        EXPECT_TRUE(parser.conditions().contains(bc::CLOAD, original));
        EXPECT_DOUBLE_EQ(original->values_start_[0], Precision(100));
        EXPECT_DOUBLE_EQ(original->values_[0], magnitude);
    }

    // NEW affects only still-inherited input definitions.
    parser.clear_conditions(bc::CLOAD);
    EXPECT_DOUBLE_EQ(original->values_[0], Precision(300));
}

TEST(Reader_Conditions, ModRetiresOldAndAddsNewSourceWhenOrientationChanges) {
    io::reader::Parser parser;
    auto old_region = std::make_shared<model::NodeRegion>("OLD");
    auto new_region = std::make_shared<model::NodeRegion>("NEW");
    old_region->add(0);
    new_region->add(1);

    auto original = std::make_shared<bc::CLoad>();
    original->region_ = old_region;
    original->values_[0] = Precision(100);
    parser.add_condition(bc::CLOAD, original, "NSET:SOURCE");

    parser.begin_loadcase(std::make_unique<loadcase::LinearStatic>());
    original->values_start_[0] = Precision(100);

    auto incoming = std::make_shared<bc::CLoad>();
    incoming->region_ = new_region;
    incoming->values_[0] = Precision(200);
    parser.modify_conditions(bc::CLOAD, "NSET:SOURCE", {incoming});

    EXPECT_TRUE(parser.conditions().contains(bc::CLOAD, original));
    EXPECT_TRUE(parser.conditions().contains(bc::CLOAD, incoming));
    EXPECT_DOUBLE_EQ(original->values_start_[0], Precision(100));
    EXPECT_DOUBLE_EQ(original->values_[0], Precision(0));
    EXPECT_DOUBLE_EQ(incoming->values_start_[0], Precision(0));
    EXPECT_DOUBLE_EQ(incoming->values_[0], Precision(200));

    parser.clear_conditions(bc::CLOAD);
    EXPECT_DOUBLE_EQ(incoming->values_[0], Precision(200));
}

TEST(Reader_Conditions, ModMatchesOriginalSourceNotOverlappingRegion) {
    io::reader::Parser parser;
    auto region = std::make_shared<model::NodeRegion>("SHARED");
    region->add(0);

    auto first = std::make_shared<bc::CLoad>();
    first->region_ = region;
    first->values_[0] = Precision(100);
    parser.add_condition(bc::CLOAD, first, "NSET:LEFT");

    auto second = std::make_shared<bc::CLoad>();
    second->region_ = region;
    second->values_[0] = Precision(200);
    parser.add_condition(bc::CLOAD, second, "NSET:RIGHT");

    parser.begin_loadcase(std::make_unique<loadcase::LinearStatic>());

    auto update = std::make_shared<bc::CLoad>();
    update->region_ = region;
    update->values_[0] = Precision(300);
    parser.modify_conditions(bc::CLOAD, "NSET:LEFT", {update});

    EXPECT_EQ(parser.conditions().get(bc::CLOAD).size(), 2u);
    EXPECT_DOUBLE_EQ(first->values_[0], Precision(300));
    EXPECT_DOUBLE_EQ(second->values_[0], Precision(200));

    parser.clear_conditions(bc::CLOAD);
    EXPECT_DOUBLE_EQ(first->values_[0], Precision(300));
    EXPECT_DOUBLE_EQ(second->values_[0], Precision(0));
}

TEST(Reader_Conditions, ModRetiresAmplitudeFromItsEffectiveInitialValue) {
    io::reader::Parser parser;
    auto region = std::make_shared<model::NodeRegion>("SOURCE");
    region->add(0);

    auto amp = std::make_shared<bc::Amplitude>("HISTORY");
    amp->add_sample(Precision(0), Precision(0.25));
    amp->add_sample(Precision(1), Precision(1));

    auto old_load = std::make_shared<bc::CLoad>();
    old_load->region_ = region;
    old_load->values_[0] = Precision(100);
    old_load->amplitude_ = amp;
    parser.add_condition(bc::CLOAD, old_load, "NSET:SOURCE");
    parser.begin_loadcase(std::make_unique<loadcase::LinearStatic>());
    old_load->values_start_[0] = Precision(100);

    auto new_load = std::make_shared<bc::CLoad>();
    new_load->region_ = region;
    new_load->values_[0] = Precision(200);
    parser.modify_conditions(bc::CLOAD, "NSET:SOURCE", {new_load});

    EXPECT_TRUE(parser.conditions().contains(bc::CLOAD, old_load));
    EXPECT_TRUE(parser.conditions().contains(bc::CLOAD, new_load));
    EXPECT_EQ(old_load->amplitude_, nullptr);
    EXPECT_DOUBLE_EQ(old_load->values_start_[0], Precision(25));
    EXPECT_DOUBLE_EQ(old_load->values_[0], Precision(0));
    EXPECT_DOUBLE_EQ(new_load->values_start_[0], Precision(0));
}

TEST(Reader_Conditions, NewPreservesTemperatureDefinedInCurrentStep) {
    io::reader::Parser parser;

    auto previous_region = std::make_shared<model::NodeRegion>("PREVIOUS");
    auto current_region  = std::make_shared<model::NodeRegion>("CURRENT");
    auto inherited = std::make_shared<bc::Temperature>(previous_region, Precision(300));
    auto current   = std::make_shared<bc::Temperature>(current_region, Precision(350));

    parser.add_condition(bc::TEMPERATURE, inherited);
    parser.begin_loadcase(std::make_unique<loadcase::LinearStatic>());
    parser.add_condition(bc::TEMPERATURE, current);

    parser.clear_conditions(bc::TEMPERATURE);

    // Prescribed-temperature constraints are removed rather than ramped down.
    EXPECT_FALSE(parser.conditions().contains(bc::TEMPERATURE, inherited));
    EXPECT_TRUE(parser.conditions().contains(bc::TEMPERATURE, current));
}

TEST(Reader_Parser, ParsesRbmCommand) {
    const std::string input_path = "tests/TMP_RBM.INP";
    const std::string output_path = "tests/TMP_RBM.RES";

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);

    {
        std::ofstream os(input_path);
        ASSERT_TRUE(os.is_open());
        os << "*NODE\n";
        os << "1, 0.0, 0.0, 0.0\n";
        os << "2, 1.0, 0.0, 0.0\n";
        os << "3, 0.0, 1.0, 0.0\n";
        os << "4, 0.0, 0.0, 1.0\n";
        os << "*ELEMENT, TYPE=C3D4, ELSET=SOLID\n";
        os << "1, 1, 2, 3, 4\n";
        os << "*MATERIAL, NAME=MAT\n*ELASTIC\n1000., 0.3\n*DENSITY\n1.\n";
        os << "*SOLID SECTION, ELSET=SOLID, MATERIAL=MAT\n";
        os << "*RBM, ELSET=SOLID\n";
    }

    io::reader::Parser parser;
    ASSERT_NO_THROW(parser.run(input_path, output_path));
    ASSERT_EQ(parser.model()._data->rbms.size(), 1u);

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);
}

TEST(Reader_Parser, ParsesOrthotropicEngineeringConstants) {
    const std::string input_path = "tests/TMP_ORTHO.INP";
    const std::string output_path = "tests/TMP_ORTHO.RES";

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);

    {
        std::ofstream os(input_path);
        ASSERT_TRUE(os.is_open());
        os << "*NODE\n1, 0., 0., 0.\n";
        os << "*MATERIAL, NAME=ORTHO\n";
        os << "*ELASTIC, TYPE=ENGINEERINGCONSTANTS\n";
        os << "100.0, 200.0, 300.0, 0.12, 0.13, 0.23, 12.0, 13.0, 23.0\n";
    }

    io::reader::Parser parser;
    ASSERT_NO_THROW(parser.run(input_path, output_path));

    auto mat = parser.model()._data->materials.get("ORTHO");
    ASSERT_NE(mat, nullptr);
    auto* ortho = mat->elasticity()->as<material::OrthotropicElasticity>();
    ASSERT_NE(ortho, nullptr);

    EXPECT_DOUBLE_EQ(ortho->E1, 100.0);
    EXPECT_DOUBLE_EQ(ortho->E2, 200.0);
    EXPECT_DOUBLE_EQ(ortho->E3, 300.0);
    EXPECT_DOUBLE_EQ(ortho->nu12, 0.12);
    EXPECT_DOUBLE_EQ(ortho->nu13, 0.13);
    EXPECT_DOUBLE_EQ(ortho->nu23, 0.23);
    EXPECT_DOUBLE_EQ(ortho->G12, 12.0);
    EXPECT_DOUBLE_EQ(ortho->G13, 13.0);
    EXPECT_DOUBLE_EQ(ortho->G23, 23.0);

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);
}

TEST(Materials_Orthotropic, TransverseShellShearUsesXzThenYz) {
    material::OrthotropicElasticity ortho(
        100.0, 200.0, 300.0,
        0.12, 0.13, 0.23,
        12.0, 13.0, 23.0
    );

    ShellMaterialStrainGreenLagrange strain;
    ShellMaterialStressPK2             stress;
    Mat5                               tangent;
    Precision old_state = Precision(0);
    Precision new_state = Precision(0);
    ortho.evaluate(strain, &old_state, &new_state, stress, &tangent);

    const Mat2 shear = tangent.template block<2, 2>(3, 3);

    EXPECT_NEAR(shear(0, 0), 13.0, 1e-12);
    EXPECT_NEAR(shear(1, 1), 23.0, 1e-12);
    EXPECT_NEAR(shear(0, 1), 0.0, 1e-12);
    EXPECT_NEAR(shear(1, 0), 0.0, 1e-12);
}
