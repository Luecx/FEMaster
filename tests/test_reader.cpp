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
#include "../src/material/orthotropic_elasticity.h"
#include "../src/material/strain/shell_material_strain_linearized.h"
#include "../src/material/stress/shell_material_stress_cauchy.h"
#include "../src/model/model.h"

#include <filesystem>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

using namespace fem;

TEST(Reader_Parser, RegistersRbmCommand) {
    io::reader::Parser parser;
    EXPECT_NE(parser.registry().find("RBM"), nullptr);
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

TEST(Reader_Parser, MapsOrthotropicElasticNu13ToNu31) {
    const std::string input_path = "tests/TMP_ORTHO.INP";
    const std::string output_path = "tests/TMP_ORTHO.RES";

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);

    {
        std::ofstream os(input_path);
        ASSERT_TRUE(os.is_open());
        os << "*MATERIAL, NAME=ORTHO\n";
        os << "*ELASTIC, TYPE=ORTHOTROPIC\n";
        os << "100.0, 200.0, 300.0, 23.0, 13.0, 12.0, 0.23, 0.13, 0.12\n";
    }

    io::reader::Parser parser;
    ASSERT_NO_THROW(parser.run(input_path, output_path));

    auto mat = parser.model()._data->materials.get("ORTHO");
    ASSERT_NE(mat, nullptr);
    auto* ortho = mat->elasticity()->as<material::OrthotropicElasticity>();
    ASSERT_NE(ortho, nullptr);

    EXPECT_DOUBLE_EQ(ortho->Gyz, 23.0);
    EXPECT_DOUBLE_EQ(ortho->Gzx, 13.0);
    EXPECT_DOUBLE_EQ(ortho->Gxy, 12.0);
    EXPECT_DOUBLE_EQ(ortho->vyz, 0.23);
    EXPECT_DOUBLE_EQ(ortho->vzx, 0.13 * 300.0 / 100.0);
    EXPECT_DOUBLE_EQ(ortho->vxy, 0.12);

    std::filesystem::remove(input_path);
    std::filesystem::remove(output_path);
}

TEST(Materials_Orthotropic, TransverseShellShearUsesXzThenYz) {
    material::OrthotropicElasticity ortho(
        100.0, 200.0, 300.0,
        23.0, 13.0, 12.0,
        0.23, 0.13, 0.12
    );

    ShellMaterialStrainLinearized strain;
    ShellMaterialStressCauchy     stress;
    Mat5                          tangent;
    ortho.evaluate(strain, nullptr, stress, tangent);

    const Mat2 shear = tangent.template block<2, 2>(3, 3);

    EXPECT_NEAR(shear(0, 0), 13.0, 1e-12);
    EXPECT_NEAR(shear(1, 1), 23.0, 1e-12);
    EXPECT_NEAR(shear(0, 1), 0.0, 1e-12);
    EXPECT_NEAR(shear(1, 0), 0.0, 1e-12);
}
