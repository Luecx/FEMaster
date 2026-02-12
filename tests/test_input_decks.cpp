#include "../src/input_decks/parser.h"
#include "../src/model/model.h"

#include <filesystem>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

using namespace fem;

TEST(InputDecks_Parser, RegistersRbmCommand) {
    input_decks::Parser parser;
    EXPECT_NE(parser.registry().find("RBM"), nullptr);
}

TEST(InputDecks_Parser, ParsesRbmCommand) {
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

    input_decks::Parser parser;
    ASSERT_NO_THROW(parser.run(input_path, output_path));
    ASSERT_EQ(parser.model()._data->rbms.size(), 1u);

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);
}
