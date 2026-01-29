#include "../src/dsl/line.h"
#include "../src/dsl/file.h"
#include "../src/dsl/keys.h"
#include "../src/dsl/registry.h"

#include <gtest/gtest.h>
#include <fstream>
#include <string>
#include <filesystem>

using namespace fem;

// 5) Line classification
TEST(DSL_Line, Classification) {
    dsl::Line l;
    l = std::string("// comment");
    EXPECT_EQ(l.type(), dsl::COMMENT);
    l = std::string("   # also comment");
    EXPECT_EQ(l.type(), dsl::COMMENT);
    l = std::string("   ");
    EXPECT_EQ(l.type(), dsl::EMPTY_LINE);
    l = std::string("*NODE, NSET=ALL");
    EXPECT_EQ(l.type(), dsl::KEYWORD_LINE);
    EXPECT_EQ(l.command(), std::string("NODE"));
    l = std::string("1, 2, 3");
    EXPECT_EQ(l.type(), dsl::DATA_LINE);
}

// 6) Keyword parsing and 7) Data tokenization
TEST(DSL_Line, KeywordAndDataParsing) {
    dsl::Line l;
    l = std::string("*ELASTIC, TYPE=iso, E=210e9, FLAG");
    EXPECT_EQ(l.type(), dsl::KEYWORD_LINE);
    EXPECT_EQ(l.command(), std::string("ELASTIC"));
    {
        auto keys = dsl::Keys::from_keyword_line(l);
        EXPECT_TRUE(keys.has("TYPE"));
        EXPECT_TRUE(keys.has("E"));
        EXPECT_TRUE(keys.has("FLAG"));
        EXPECT_EQ(keys.raw("TYPE"), std::string("ISO")); // uppercased
        EXPECT_EQ(keys.raw("FLAG"), std::string(""));    // empty value for flag
    }

    l = std::string(" foo,Bar , baz ");
    ASSERT_EQ(l.type(), dsl::DATA_LINE);
    const auto& vals = l.values();
    ASSERT_EQ(vals.size(), 3u);
    EXPECT_EQ(vals[0], std::string("FOO"));
    EXPECT_EQ(vals[1], std::string("BAR"));
    EXPECT_EQ(vals[2], std::string("BAZ"));
}

// 8) File includes: *INCLUDE,SRC=...
TEST(DSL_File, Include) {
    // Note: The DSL uppercases keyword lines, including file paths in INCLUDE.
    // Use uppercase paths on disk to keep lookups working on case-sensitive FS.
    const std::string inc = "TESTS/TMP_INC.DECK";
    const std::string mainf = "TESTS/TMP_MAIN.DECK";
    std::filesystem::create_directories("TESTS");
    {
        std::ofstream os(inc);
        os << "1,2,3\n";
        os << "*FOO,BAR\n";
    }
    {
        std::ofstream os(mainf);
        os << "*INCLUDE,SRC=" << inc << "\n";
        os << "4,5,6\n";
    }

    dsl::File f(mainf);
    auto& l1 = f.next_line();
    ASSERT_EQ(l1.type(), dsl::DATA_LINE);
    EXPECT_EQ(l1.values().size(), 3u);
    auto& l2 = f.next_line();
    ASSERT_EQ(l2.type(), dsl::KEYWORD_LINE);
    EXPECT_EQ(l2.command(), std::string("FOO"));
    auto& l3 = f.next_line();
    ASSERT_EQ(l3.type(), dsl::DATA_LINE);
    EXPECT_EQ(l3.values().size(), 3u);
}

// 9–10) Keys normalization and boolean semantics; 11–12) Registry basics/printing
TEST(DSL_Keys_Registry, KeysAndRegistry) {
    dsl::Line l; l = std::string("*CMD, A=1, B=YES, C=off, FLAG");
    auto K = dsl::Keys::from_keyword_line(l);

    // Build spec with alias + defaults + allowed
    dsl::KeywordSpec spec = dsl::KeywordSpec::make()
        .key("ALPHA").alternative("A").allowed({"1","2"})
        .key("BETA").alternative("B").allowed({"YES","NO"})
        .flag("FLAG")
        .key("GAMMA").optional("42");

    // Apply spec (normalizes aliases, validates allowed, injects defaults)
    K.apply_spec(spec, "CMD");
    EXPECT_TRUE(K.has("ALPHA"));
    EXPECT_TRUE(K.has("BETA"));
    EXPECT_TRUE(K.has("GAMMA"));
    EXPECT_TRUE(K.has("FLAG"));
    EXPECT_EQ(K.get<std::string>("GAMMA"), std::string("42"));

    // Boolean semantics
    EXPECT_TRUE(K.get<bool>("FLAG"));
    EXPECT_TRUE(K.get<bool>("BETA"));         // YES → true
    // Introduce an explicit false token
    dsl::Line l2; l2 = std::string("*CMD, X=OFF, Y=0, Z=N");
    auto K2 = dsl::Keys::from_keyword_line(l2);
    EXPECT_FALSE(K2.get<bool>("X"));
    EXPECT_FALSE(K2.get<bool>("Y"));
    EXPECT_FALSE(K2.get<bool>("Z"));

    // Registry create/find and print helpers sanity
    dsl::Registry reg;
    reg.command("ELASTIC", [](dsl::Command&){ /* no-op */ });
    ASSERT_NE(reg.find("ELASTIC"), nullptr);
    // Smoke test for printing (no throw)
    reg.print_help();
}
