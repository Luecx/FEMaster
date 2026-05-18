#include "../src/data/field.h"
#include "../src/reader/writer.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

using namespace fem;

namespace {

std::string read_text(const std::string& path) {
    std::ifstream input(path);
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

} // namespace

TEST(Reader_Writer, WritesFieldTypeForModelField) {
    const std::string output_path = "tests/TMP_WRITER_FIELD.RES";
    std::filesystem::remove(output_path);

    model::Field field("U", model::FieldDomain::NODE, 2, 3);
    field.set_zero();

    {
        reader::Writer writer(output_path);
        writer.write_field(field, "DISPLACEMENT");
    }

    const std::string text = read_text(output_path);
    EXPECT_NE(text.find("FIELD, NAME=DISPLACEMENT, TYPE=NODE, COLS=3, ROWS=2"), std::string::npos);

    std::filesystem::remove(output_path);
}

TEST(Reader_Writer, WritesInferredTypeForIndexedMatrixField) {
    const std::string output_path = "tests/TMP_WRITER_MATRIX.RES";
    std::filesystem::remove(output_path);

    DynamicMatrix matrix(1, 5);
    matrix << 1.0, 2.0, 10.0, 20.0, 30.0;

    {
        reader::Writer writer(output_path);
        writer.write_eigen_matrix(matrix, "LOCAL_SECTION_FORCES", 2);
    }

    const std::string text = read_text(output_path);
    EXPECT_NE(text.find("FIELD, NAME=LOCAL_SECTION_FORCES, TYPE=ELEMENT_NODAL, INDEX_COLS=2, VALUE_COLS=3, ROWS=1"),
              std::string::npos);

    std::filesystem::remove(output_path);
}
