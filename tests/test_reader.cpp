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
#include "../src/io/reader/parser_abq.h"
#include "../src/io/dsl/deck_parser.h"
#include "../src/io/dsl/file.h"
#include "../src/bc/neumann/load_c.h"
#include "../src/loadcase/linear_static.h"
#include "../src/bc/neumann/load_inertial.h"
#include "../src/material/orthotropic_elasticity.h"
#include "../src/material/strain/shell_material_strain_linearized.h"
#include "../src/material/stress/shell_material_stress_cauchy.h"
#include "../src/model/model.h"

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <stdexcept>
#include <initializer_list>

#include <gtest/gtest.h>

using namespace fem;

TEST(Reader_Parser, RegistersRbmCommand) {
    io::reader::Parser parser;
    EXPECT_NE(parser.registry().find("RBM"), nullptr);
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
    ASSERT_EQ(collector->entries().size(), 1u);

    auto load = std::dynamic_pointer_cast<bc::InertialLoad>(collector->entries().front());
    ASSERT_NE(load, nullptr);
    EXPECT_EQ(load->region_, parser.model()._data->elem_sets.get("EALL"));
    EXPECT_DOUBLE_EQ(load->center_acc_(2), 9.81);

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);
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
        os << "*RBM, NSET=NALL\n";
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

    ShellMaterialStrainLinearized strain;
    ShellMaterialStressCauchy     stress;
    Mat5                          tangent;
    Precision old_state = Precision(0);
    Precision new_state = Precision(0);
    ortho.evaluate(strain, &old_state, &new_state, stress, tangent);

    const Mat2 shear = tangent.template block<2, 2>(3, 3);

    EXPECT_NEAR(shear(0, 0), 13.0, 1e-12);
    EXPECT_NEAR(shear(1, 1), 23.0, 1e-12);
    EXPECT_NEAR(shear(0, 1), 0.0, 1e-12);
    EXPECT_NEAR(shear(1, 0), 0.0, 1e-12);
}

namespace {
struct TempCloadDeck {
    std::string file;
    TempCloadDeck(const std::string& name, const std::string& contents)
        : file("tests/" + name) {
        std::ofstream stream(file);
        if (!stream) throw std::runtime_error("Cannot create temporary CLOAD deck");
        stream << contents;
    }
    ~TempCloadDeck() { std::filesystem::remove(file); }
};
} // namespace

TEST(Reader_CLoad, NativeAcceptsMixedAbaqusAndVectorRowsWithoutNewFiles) {
    TempCloadDeck deck("TMP_SHARED_CLOAD_NATIVE.inp",
        "*NODE\n1, 0., 0., 0.\n2, 1., 0., 0.\n"
        "*NSET, NAME=ENDS\n1, 2\n"
        "*CLOAD, LOAD_COLLECTOR=FORCES\n"
        "ENDS, 1, 100.\nENDS, 0., 5., -10.\nENDS, 6, 7.\n");
    io::reader::Parser parser;
    ASSERT_NO_THROW(parser.run(deck.file, "tests/TMP_SHARED_CLOAD_NATIVE"));
    const auto collector = parser.model()._data->load_cols.get("FORCES");
    ASSERT_NE(collector, nullptr);
    ASSERT_EQ(collector->entries().size(), 3u);
    const auto a = std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[0]);
    const auto b = std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[1]);
    const auto c = std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[2]);
    ASSERT_NE(a, nullptr);
    ASSERT_NE(b, nullptr);
    ASSERT_NE(c, nullptr);
    EXPECT_EQ(a->region_->size(), 2u);
    EXPECT_EQ(b->region_->size(), 2u);
    EXPECT_EQ(c->region_->size(), 2u);
    EXPECT_DOUBLE_EQ(a->values_[0], 100.);
    EXPECT_DOUBLE_EQ(b->values_[1], 5.);
    EXPECT_DOUBLE_EQ(b->values_[2], -10.);
    EXPECT_DOUBLE_EQ(c->values_[5], 7.);
}

TEST(Reader_CLoad, AbaqusReaderUsesSameRegisteredCloadForNamedModelLoads) {
    TempCloadDeck deck("TMP_SHARED_CLOAD_ABQ.inp",
        "*NODE\n1, 0., 0., 0.\n"
        "*CLOAD, LOAD_COLLECTOR=FORCES\n"
        "1, 1, 12.\n1, 1., 2., 3., 4., 5., 6.\n");
    io::reader::ParserAbq parser;
    ASSERT_NO_THROW(parser.run(deck.file, "tests/TMP_SHARED_CLOAD_ABQ"));
    const auto collector = parser.model()._data->load_cols.get("FORCES");
    ASSERT_NE(collector, nullptr);
    ASSERT_EQ(collector->entries().size(), 2u);
    EXPECT_DOUBLE_EQ(std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[0])->values_[0], 12.);
    EXPECT_DOUBLE_EQ(std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[1])->values_[5], 6.);
}

TEST(Reader_CLoad, UnnamedLoadsOutsideAnalysesAreRejectedByBothReaders) {
    TempCloadDeck deck("TMP_SHARED_CLOAD_UNNAMED.inp",
        "*NODE\n1, 0., 0., 0.\n*CLOAD\n1, 1, 12.\n");
    io::reader::Parser native;
    io::reader::ParserAbq abq;
    EXPECT_THROW(native.run(deck.file, "tests/TMP_SHARED_CLOAD_UNNAMED_NATIVE"), std::exception);
    EXPECT_THROW(abq.run(deck.file, "tests/TMP_SHARED_CLOAD_UNNAMED_ABQ"), std::exception);
}

TEST(Reader_CLoad, NativeLoadcaseAndStepInlineCollectors) {
    for (const bool use_step : {false, true}) {
        TempCloadDeck deck(use_step ? "TMP_SHARED_CLOAD_STEP.inp" : "TMP_SHARED_CLOAD_CASE.inp",
            use_step
                ? "*STEP\n*STATIC\n*CLOAD\n1, 1, 25.\n1, 0., 0., -5.\n*END STEP\n"
                : "*LOADCASE, TYPE=LINEARSTATIC\n*CLOAD\n1, 1, 25.\n1, 0., 0., -5.\n*END\n");
        io::reader::Parser parser;
        parser.model().set_node(1, 0., 0., 0.);
        parser.model().compile();
        io::dsl::File file(deck.file);
        const auto parsed = io::dsl::DeckParser(parser.registry()).parse(file);
        const auto scopes = parsed.root().children(use_step ? "STEP" : "LOADCASE");
        ASSERT_EQ(scopes.size(), 1u);
        scopes.front()->enter();
        if (use_step) scopes.front()->execute_children("STATIC");
        scopes.front()->execute_children("CLOAD");

        const auto* active = dynamic_cast<loadcase::LinearStatic*>(parser.active_loadcase());
        ASSERT_NE(active, nullptr);
        ASSERT_EQ(active->loads.size(), 1u);
        const auto collector = parser.model()._data->load_cols.get(active->loads.front());
        ASSERT_NE(collector, nullptr);
        ASSERT_EQ(collector->entries().size(), 2u);
        EXPECT_DOUBLE_EQ(std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[0])->values_[0], 25.);
        EXPECT_DOUBLE_EQ(std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[1])->values_[2], -5.);
    }
}

TEST(Reader_CLoad, AbaqusStepUsesInlineCollectorAndSharedSyntax) {
    TempCloadDeck deck("TMP_SHARED_CLOAD_ABQ_STEP.inp",
        "*STEP\n*STATIC\n*CLOAD\n1, 1, 25.\n1, 0., 0., -5.\n*END STEP\n");
    io::reader::ParserAbq parser;
    parser.model().set_node(1, 0., 0., 0.);
    parser.model().compile();
    io::dsl::File file(deck.file);
    const auto parsed = io::dsl::DeckParser(parser.registry()).parse(file);
    const auto scopes = parsed.root().children("STEP");
    ASSERT_EQ(scopes.size(), 1u);
    scopes.front()->enter();
    scopes.front()->execute_children("STATIC");
    scopes.front()->execute_children("CLOAD");

    const auto* active = dynamic_cast<loadcase::LinearStatic*>(parser.active_loadcase());
    ASSERT_NE(active, nullptr);
    ASSERT_EQ(active->loads.size(), 1u);
    const auto collector = parser.model()._data->load_cols.get(active->loads.front());
    ASSERT_NE(collector, nullptr);
    ASSERT_EQ(collector->entries().size(), 2u);
    EXPECT_DOUBLE_EQ(std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[0])->values_[0], 25.);
    EXPECT_DOUBLE_EQ(std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[1])->values_[2], -5.);
}

TEST(Reader_CLoad, ShortVectorRowsCannotBeMistakenForAbaqusDofRows) {
    TempCloadDeck deck("TMP_SHARED_CLOAD_BAD_VECTOR.inp",
        "*NODE\n1, 0., 0., 0.\n"
        "*CLOAD, LOAD_COLLECTOR=FORCES\n"
        "1, 1.0, 100.\n");
    io::reader::Parser parser;
    EXPECT_THROW(parser.run(deck.file, "tests/TMP_SHARED_CLOAD_BAD_VECTOR"), std::exception);
}


TEST(Reader_CLoad, SequentialNativeLoadcasesAndStepKeepIndependentLoads) {
    // A complete three-bar truss: this exercises solver execution and prevents
    // one loadcase from silently reusing the other's implicit load collector.
    TempCloadDeck deck("TMP_SHARED_CLOAD_TWO_CASES.inp",
        "*NODE\n1, 0., 0., 0.\n2, 1000., 0., 0.\n3, 500., 300., 0.\n"
        "*ELEMENT, TYPE=T3, ELSET=TRUSS\n1, 1, 3\n2, 2, 3\n3, 1, 2\n"
        "*NSET, NAME=LEFT\n1\n*NSET, NAME=RIGHT\n2\n*NSET, NAME=TOP\n3\n"
        "*MATERIAL, NAME=STEEL\n*ELASTIC, TYPE=ISOTROPIC\n210000., 0.3\n"
        "*TRUSSSECTION, ELSET=TRUSS, MATERIAL=STEEL\n100.\n"
        "*SUPPORT, SUPPORT_COLLECTOR=BC\n"
        "LEFT, 0., 0., 0., NAN, NAN, NAN\n"
        "RIGHT, NAN, 0., 0., NAN, NAN, NAN\n"
        "TOP, NAN, NAN, 0., NAN, NAN, NAN\n"
        "*LOADCASE, TYPE=LINEARSTATIC\n*SUPPORTS\nBC\n"
        "*CLOAD\nTOP, 1, 100.\n*END\n"
        "*LOADCASE, TYPE=LINEARSTATIC\n*SUPPORTS\nBC\n"
        "*CLOAD\nTOP, 1, 200.\n*END\n"
        "*STEP\n*STATIC\n*SUPPORTS\nBC\n*CLOAD\nTOP, 1, 300.\n*END STEP\n");
    io::reader::Parser parser;
    ASSERT_NO_THROW(parser.run(deck.file, "tests/TMP_SHARED_CLOAD_TWO_CASES"));
    const ID top = parser.model().compiled_node_id("3");
    const auto first = parser.model().build_load_matrix({"__FEMASTER_INLINE_CLOAD_1"});
    const auto second = parser.model().build_load_matrix({"__FEMASTER_INLINE_CLOAD_2"});
    EXPECT_DOUBLE_EQ(first(top, 0), 100.);
    EXPECT_DOUBLE_EQ(second(top, 0), 200.);
    const auto step_loads = parser.model().build_load_matrix({"__ABQ_STEP_LOADS"});
    EXPECT_DOUBLE_EQ(step_loads(top, 0), 300.);
    EXPECT_DOUBLE_EQ(first(top, 1), 0.);
    EXPECT_DOUBLE_EQ(second(top, 1), 0.);
}

TEST(Reader_CLoad, NamedCollectorReferencedInsideLoadcaseIsAppliedOnce) {
    TempCloadDeck deck("TMP_SHARED_CLOAD_COLLECTOR_DEDUP.inp",
        "*LOADCASE, TYPE=LINEARSTATIC\n"
        "*CLOAD, LOAD_COLLECTOR=SHARED\n1, 1, 10.\n"
        "*LOADS\nSHARED\n"
        "*CLOAD, LOAD_COLLECTOR=SHARED\n1, 2, 20.\n"
        "*CLOAD\n1, 0., 0., 30.\n*END\n");
    io::reader::Parser parser;
    parser.model().set_node(1, 0., 0., 0.);
    parser.model().compile();
    io::dsl::File file(deck.file);
    const auto parsed = io::dsl::DeckParser(parser.registry()).parse(file);
    const auto scopes = parsed.root().children("LOADCASE");
    ASSERT_EQ(scopes.size(), 1u);
    scopes.front()->enter();
    scopes.front()->execute_children("CLOAD");
    scopes.front()->execute_children("LOADS");

    const auto* active = dynamic_cast<loadcase::LinearStatic*>(parser.active_loadcase());
    ASSERT_NE(active, nullptr);
    ASSERT_EQ(active->loads.size(), 2u);
    const ID node = parser.model().compiled_node_id("1");
    const auto forces = parser.model().build_load_matrix(active->loads);
    EXPECT_DOUBLE_EQ(forces(node, 0), 10.);
    EXPECT_DOUBLE_EQ(forces(node, 1), 20.);
    EXPECT_DOUBLE_EQ(forces(node, 2), 30.);
}
