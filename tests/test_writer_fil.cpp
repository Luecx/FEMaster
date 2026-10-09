/**
 * @file test_writer_fil.cpp
 * @brief ASCII .fil record framing and first-frame results smoke tests.
 */
#include "../src/io/writer/writer_fil.h"
#include "../src/model/model.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include <gtest/gtest.h>

using namespace fem;

namespace {
std::string fil_contents(const std::string& filename) {
    std::ifstream stream(filename, std::ios::binary);
    return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
}
}

TEST(FilWriter, WritesModelAndTwoFramesWithEightyColumnRecords) {
    const std::string path = "tests/TMP_FIL_WRITER.fil";
    std::filesystem::remove(path);
    model::Model model;
    model.set_node(0, 0., 0., 0.);
    model.set_node(1, 1., 2., 3.);
    model.compile();

    model::Field displacement("U", model::FieldDomain::NODE, 2, 3);
    displacement.set_zero();
    displacement(1, 0) = 0.25;

    {
        io::writer::FilWriter writer(path);
        writer.write_model_data(*model._data);
        writer.add_loadcase(2, io::writer::WriterStepType::Static);
        writer.begin_frame(0.25);
        writer.write_field(displacement, "DISPLACEMENT", model._data.get(), 0.25);
        writer.end_frame();
        writer.begin_frame(0.5);
        writer.write_field(displacement, "DISPLACEMENT", model._data.get(), 0.5);
        writer.end_frame();
    }

    const std::string text = fil_contents(path);
    EXPECT_NE(text.find("I 41921"), std::string::npos);
    EXPECT_NE(text.find("I 41901"), std::string::npos);
    EXPECT_NE(text.find("I 42000"), std::string::npos);
    EXPECT_NE(text.find("I 3101"), std::string::npos);
    EXPECT_NE(text.find("I 42001"), std::string::npos);
    EXPECT_NE(text.find("D 2.500000000000000D-01"), std::string::npos);

    std::size_t begin = 0;
    int lines = 0;
    while (begin < text.size()) {
        const auto end = text.find('\n', begin);
        ASSERT_NE(end, std::string::npos);
        EXPECT_EQ(end - begin, 80u);
        ++lines;
        begin = end + 1;
    }
    EXPECT_GT(lines, 2);
    const std::string terminator(80, ' ');
    EXPECT_NE(text.find(terminator + "\n"), std::string::npos);
    std::filesystem::remove(path);
}
