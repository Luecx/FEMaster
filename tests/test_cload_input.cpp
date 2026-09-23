/**
 * @file test_cload_input.cpp
 * @brief Regression coverage for shared CLOAD syntax and collector scoping.
 */
#include "../src/io/reader/commands/cload_common.h"
#include "../src/io/reader/parser.h"
#include "../src/io/reader/parser_abq.h"
#include "../src/io/dsl/deck_parser.h"
#include "../src/io/dsl/file.h"
#include "../src/loadcase/linear_static.h"
#include "../src/model/model.h"

#include <array>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <memory>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

using namespace fem;

namespace {
std::array<std::string, 6> row(std::initializer_list<std::string> given) {
    std::array<std::string, 6> values;
    values.fill(io::reader::commands::cload_common::missing_token);
    std::size_t index = 0;
    for (const auto& value : given) values.at(index++) = value;
    return values;
}

struct TemporaryDeck {
    std::string file;
    explicit TemporaryDeck(const std::string& filename, const std::string& deck)
        : file("tests/" + filename) {
        std::ofstream stream(file);
        if (!stream) throw std::runtime_error("Cannot open temporary input deck");
        stream << deck;
    }
    ~TemporaryDeck() { std::filesystem::remove(file); }
};
} // namespace

TEST(CLoad_Input, DOFAndVectorRowsHaveUnambiguousLengths) {
    const auto dof = io::reader::commands::cload_common::parse(row({"3", "-500."}));
    EXPECT_DOUBLE_EQ(dof[2], -500.);
    EXPECT_DOUBLE_EQ(dof[0], 0.);

    const auto vec = io::reader::commands::cload_common::parse(row({"3", "-500.", "2."}));
    EXPECT_DOUBLE_EQ(vec[0], 3.);
    EXPECT_DOUBLE_EQ(vec[1], -500.);
    EXPECT_DOUBLE_EQ(vec[2], 2.);

    const auto moment = io::reader::commands::cload_common::parse(row({"6", "40."}));
    EXPECT_DOUBLE_EQ(moment[5], 40.);
    EXPECT_DOUBLE_EQ(moment[0], 0.);
    EXPECT_THROW(io::reader::commands::cload_common::parse(row({"7", "40."})), std::exception);
    EXPECT_THROW(io::reader::commands::cload_common::parse(row({"1.0", "40."})), std::exception);
    EXPECT_THROW(io::reader::commands::cload_common::parse(row({"100."})), std::exception);
    EXPECT_THROW(io::reader::commands::cload_common::parse(row({"1", "2", "INVALID"})), std::exception);
}

TEST(CLoad_Input, NativeModelLevelBlockAcceptsMixedRowsAndNodeSets) {
    TemporaryDeck input("TMP_CLOAD_NATIVE.inp",
        "*NODE\n1, 0., 0., 0.\n2, 1., 0., 0.\n"
        "*NSET, NAME=ENDS\n1, 2\n"
        "*CLOAD, LOAD_COLLECTOR=FORCES\n"
        "ENDS, 1, 100.\nENDS, 0., 5., -10.\nENDS, 6, 7.\n");
    io::reader::Parser parser;
    ASSERT_NO_THROW(parser.run(input.file, "tests/TMP_CLOAD_NATIVE"));
    const auto collector = parser.model()._data->load_cols.get("FORCES");
    ASSERT_NE(collector, nullptr);
    ASSERT_EQ(collector->entries().size(), 3u);
    for (const auto& load : collector->entries()) {
        const auto concentrated = std::dynamic_pointer_cast<bc::CLoad>(load);
        ASSERT_NE(concentrated, nullptr);
        EXPECT_EQ(concentrated->region_->size(), 2u);
    }
    EXPECT_DOUBLE_EQ(std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[0])->values_[0], 100.);
    EXPECT_DOUBLE_EQ(std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[1])->values_[1], 5.);
    EXPECT_DOUBLE_EQ(std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[2])->values_[5], 7.);
}

TEST(CLoad_Input, BothReadersRejectUnnamedGlobalLoads) {
    const std::string deck = "*NODE\n1, 0., 0., 0.\n*CLOAD\n1, 1, 100.\n";
    TemporaryDeck native("TMP_CLOAD_NONAME_NATIVE.inp", deck);
    TemporaryDeck abq("TMP_CLOAD_NONAME_ABQ.inp", deck);
    io::reader::Parser native_parser;
    io::reader::ParserAbq abq_parser;
    EXPECT_THROW(native_parser.run(native.file, "tests/TMP_CLOAD_NONAME_NATIVE"), std::exception);
    EXPECT_THROW(abq_parser.run(abq.file, "tests/TMP_CLOAD_NONAME_ABQ"), std::exception);
}

TEST(CLoad_Input, AbaqusReaderAcceptsNamedGlobalVectorAndDOFRows) {
    TemporaryDeck input("TMP_CLOAD_ABQ.inp",
        "*NODE\n1, 0., 0., 0.\n"
        "*CLOAD, LOAD_COLLECTOR=GLOBAL\n"
        "1, 1, 12.\n1, 1., 2., 3., 4., 5., 6.\n");
    io::reader::ParserAbq parser;
    ASSERT_NO_THROW(parser.run(input.file, "tests/TMP_CLOAD_ABQ"));
    const auto collector = parser.model()._data->load_cols.get("GLOBAL");
    ASSERT_NE(collector, nullptr);
    ASSERT_EQ(collector->entries().size(), 2u);
    EXPECT_DOUBLE_EQ(std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[0])->values_[0], 12.);
    EXPECT_DOUBLE_EQ(std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[1])->values_[5], 6.);
}

TEST(CLoad_Input, InlineCollectorsAreAutomaticallySelectedWithoutDuplication) {
    io::reader::Parser parser;
    parser.begin_loadcase(std::make_unique<loadcase::LinearStatic>());
    parser.activate_cload_collector("", true);
    auto* lc = dynamic_cast<loadcase::LinearStatic*>(parser.active_loadcase());
    ASSERT_NE(lc, nullptr);
    ASSERT_EQ(lc->loads.size(), 1u);
    EXPECT_EQ(lc->loads.front(), "__INTERNAL_CLOAD_1");

    parser.activate_cload_collector("", true);
    EXPECT_EQ(lc->loads.size(), 1u);
    parser.activate_cload_collector("OTHER", true);
    parser.activate_cload_collector("OTHER", true);
    ASSERT_EQ(lc->loads.size(), 2u);
    EXPECT_EQ(lc->loads.back(), "OTHER");
}

TEST(CLoad_Input, BothReadersParseInlineMixedRowsUnderAnalysisScopes) {
    TemporaryDeck native("TMP_CLOAD_SCOPE_NATIVE.inp",
        "*LOADCASE, TYPE=LINEARSTATIC\n"
        "*CLOAD\nTIP, 1, 100.\nTIP, 0., 0., -200.\n"
        "*END\n");
    TemporaryDeck abq("TMP_CLOAD_SCOPE_ABQ.inp",
        "*STEP\n*STATIC\n"
        "*CLOAD\nTIP, 1, 100.\nTIP, 0., 0., -200.\n"
        "*END STEP\n");

    io::reader::Parser native_parser;
    io::reader::ParserAbq abq_parser;
    io::dsl::File native_file(native.file);
    io::dsl::File abq_file(abq.file);
    const auto native_deck = io::dsl::DeckParser(native_parser.registry()).parse(native_file);
    const auto cases = native_deck.root().children("LOADCASE");
    ASSERT_EQ(cases.size(), 1u);
    ASSERT_EQ(cases.front()->children("CLOAD").size(), 1u);

    const auto abq_deck = io::dsl::DeckParser(abq_parser.registry()).parse(abq_file);
    const auto steps = abq_deck.root().children("STEP");
    ASSERT_EQ(steps.size(), 1u);
    ASSERT_EQ(steps.front()->children("CLOAD").size(), 1u);
}

TEST(CLoad_Input, NativeInlineLoadIsMaterializedAndSelectedBeforeSolving) {
    TemporaryDeck input("TMP_CLOAD_INLINE_NATIVE.inp",
        "*LOADCASE, TYPE=LINEARSTATIC\n"
        "*CLOAD\n1, 1, 25.\n1, 0., 0., -5.\n"
        "*END\n");
    io::reader::Parser parser;
    parser.model().set_node(1, 0., 0., 0.);
    parser.model().compile();

    io::dsl::File file(input.file);
    const auto deck = io::dsl::DeckParser(parser.registry()).parse(file);
    const auto cases = deck.root().children("LOADCASE");
    ASSERT_EQ(cases.size(), 1u);
    cases.front()->enter();
    cases.front()->execute_children("CLOAD"); // Do not execute the solver in this parser test.

    auto* lc = dynamic_cast<loadcase::LinearStatic*>(parser.active_loadcase());
    ASSERT_NE(lc, nullptr);
    ASSERT_EQ(lc->loads.size(), 1u);
    const auto collector = parser.model()._data->load_cols.get(lc->loads.front());
    ASSERT_NE(collector, nullptr);
    ASSERT_EQ(collector->entries().size(), 2u);
    EXPECT_DOUBLE_EQ(std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[0])->values_[0], 25.);
    EXPECT_DOUBLE_EQ(std::dynamic_pointer_cast<bc::CLoad>(collector->entries()[1])->values_[2], -5.);
}

TEST(CLoad_Input, AbaqusInlineLoadPreservesDefaultStepCollector) {
    TemporaryDeck input("TMP_CLOAD_INLINE_ABQ.inp",
        "*STEP\n*STATIC\n"
        "*CLOAD, LOAD_COLLECTOR=SPECIAL\n1, 1, 25.\n"
        "*CLOAD\n1, 0., 0., -5.\n"
        "*END STEP\n");
    io::reader::ParserAbq parser;
    parser.model().set_node(1, 0., 0., 0.);
    parser.model().compile();

    io::dsl::File file(input.file);
    const auto deck = io::dsl::DeckParser(parser.registry()).parse(file);
    const auto steps = deck.root().children("STEP");
    ASSERT_EQ(steps.size(), 1u);
    steps.front()->enter();
    steps.front()->execute_children("STATIC");
    steps.front()->execute_children("CLOAD"); // Do not execute END STEP / solver.

    auto* lc = dynamic_cast<loadcase::LinearStatic*>(parser.active_loadcase());
    ASSERT_NE(lc, nullptr);
    ASSERT_EQ(lc->loads.size(), 2u);
    const auto named = parser.model()._data->load_cols.get("SPECIAL");
    const auto defaults = parser.model()._data->load_cols.get("__ABQ_STEP_LOADS");
    ASSERT_NE(named, nullptr);
    ASSERT_NE(defaults, nullptr);
    EXPECT_EQ(named->entries().size(), 1u);
    EXPECT_EQ(defaults->entries().size(), 1u);
    EXPECT_EQ(parser.model()._data->load_cols.get(), defaults);
}
