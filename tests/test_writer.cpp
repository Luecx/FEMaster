/**
 * @file test_writer.cpp
 * @brief Tests result-writer metadata and compiled shell-resultant output.
 */

#include "../src/data/field.h"
#include "../src/material/isotropic_elasticity.h"
#include "../src/model/model.h"
#include "../src/model/shell/s4.h"
#include "../src/io/writer/writer_frd.h"
#include "../src/io/writer/writer_res.h"
#include "../src/section/section_shell_integrated.h"

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
        io::writer::ResWriter writer(output_path);
        writer.write_field(field, "DISPLACEMENT");
    }

    const std::string text = read_text(output_path);
    EXPECT_NE(text.find("FIELD, NAME=DISPLACEMENT, TYPE=NODE, COLS=3, ROWS=2"), std::string::npos);

    std::filesystem::remove(output_path);
}

TEST(Reader_Writer, WritesInferredTypeForIndexedMatrixField) {
    const std::string output_path = "tests/TMP_WRITER_MATRIX.RES";
    std::filesystem::remove(output_path);

    model::Model model;
    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_node(2, 1.0, 1.0, 0.0);
    model.set_node(3, 0.0, 1.0, 0.0);
    model.set_element<model::S4>(0, 0, 1, 2, 3);
    model.compile();

    model::Field field("LOCAL_SECTION_FORCES", model::FieldDomain::ELEMENT_NODAL, 4, 3);
    field.set_zero();

    {
        io::writer::ResWriter writer(output_path);
        writer.write_field(field, "LOCAL_SECTION_FORCES", model._data.get());
    }

    const std::string text = read_text(output_path);
    EXPECT_NE(text.find("FIELD, NAME=LOCAL_SECTION_FORCES, TYPE=ELEMENT_NODAL, INDEX_COLS=2, VALUE_COLS=3, ROWS=4"),
              std::string::npos);

    std::filesystem::remove(output_path);
}

TEST(Reader_Writer, WritesEightShellResultantComponentsToFrd) {
    const std::string output_path = "tests/TMP_WRITER_SHELL_RESULTANTS.FRD";
    std::filesystem::remove(output_path);

    fem::model::Model model;
    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_node(2, 1.0, 1.0, 0.0);
    model.set_node(3, 0.0, 1.0, 0.0);
    model.set_element<fem::model::S4>(0, 0, 1, 2, 3);

    auto material = std::make_shared<fem::material::Material>("MAT");
    material->set_elasticity<fem::material::IsotropicElasticity>(1000.0, 0.3);
    model.add_material(material);

    model.add_section(std::make_shared<fem::IntegratedShellSection>(
        material,
        model._data->parts.get()->elem_sets.get(fem::SET_ELEM_ALL),
        0.1
    ));
    model.compile();

    fem::model::Field displacement("U", fem::model::FieldDomain::NODE, 4, 6);
    displacement.set_zero();

    fem::model::Field resultants = model.compute_shell_resultants(displacement);

    {
        fem::io::writer::FrdWriter writer(output_path);
        writer.write_model_data(*model._data);
        writer.write_field(resultants, "SHELLRESULTANTS", model._data.get());
    }

    const std::string text = read_text(output_path);
    EXPECT_NE(text.find(" -4  SHR         8    1"), std::string::npos);
    EXPECT_NE(text.find(" -5  SHR8"), std::string::npos);

    std::filesystem::remove(output_path);
}
