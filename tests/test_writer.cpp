/**
 * @file test_writer.cpp
 * @brief Tests result-writer metadata and compiled shell-resultant output.
 */

#include "../src/data/field.h"
#include "../src/material/isotropic_elasticity.h"
#include "../src/model/model.h"
#include "../src/model/shell/s4.h"
#include "../src/io/writer/writer_frd.h"
#include "../src/io/writer/writer_femr.h"
#include "../src/io/writer/writer_res.h"
#include "../src/section/section_shell_integrated.h"

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

using namespace fem;

namespace {

std::string read_text(const std::string& path) {
    std::ifstream input(path);
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

template<typename T>
T read_binary_value(const std::vector<std::uint8_t>& bytes, std::size_t offset) {
    T value{};
    std::memcpy(&value, bytes.data() + offset, sizeof(T));
    return value;
}

std::vector<std::pair<std::uint32_t, double>> read_femr_frames(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    const std::vector<std::uint8_t> bytes(
        std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
    std::vector<std::pair<std::uint32_t, double>> frames;
    std::size_t offset = 0;
    constexpr std::size_t chunk_header_size = 32;
    while (offset + chunk_header_size <= bytes.size()) {
        const std::string type(reinterpret_cast<const char*>(bytes.data() + offset), 4);
        const auto stored_size = read_binary_value<std::uint64_t>(bytes, offset + 8);
        const std::size_t payload = offset + chunk_header_size;
        if (type == "FRAM") {
            const auto frame_id = read_binary_value<std::uint32_t>(bytes, payload + 4);
            const auto frame_value = read_binary_value<double>(bytes, payload + 8);
            frames.emplace_back(frame_id, frame_value);
        }
        offset = payload + static_cast<std::size_t>(stored_size);
    }
    return frames;
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

TEST(Reader_Writer, FemrStaticFramesFollowChangingFrameValues) {
    const std::string output_path = "tests/TMP_WRITER_STATIC_FRAMES.FEMR";
    std::filesystem::remove(output_path);

    model::Field field("VALUE", model::FieldDomain::UNKNOWN, 1, 1);
    field.set_zero();

    {
        io::writer::FemrWriter writer(output_path);
        writer.add_loadcase(1, io::writer::WriterStepType::Static);
        writer.write_field(field, "DISPLACEMENT_1", nullptr, 0.25);
        writer.write_field(field, "STRESS_1", nullptr, 0.25);
        writer.write_field(field, "DISPLACEMENT_2", nullptr, 0.5);
        writer.write_field(field, "STRESS_2", nullptr, 0.5);
    }

    const auto frames = read_femr_frames(output_path);
    ASSERT_EQ(frames.size(), 2u);
    EXPECT_EQ(frames[0].first, 0u);
    EXPECT_DOUBLE_EQ(frames[0].second, 0.25);
    EXPECT_EQ(frames[1].first, 1u);
    EXPECT_DOUBLE_EQ(frames[1].second, 0.5);

    std::filesystem::remove(output_path);
}
