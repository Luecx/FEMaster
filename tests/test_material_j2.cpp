#include "../src/io/dsl/deck_parser.h"
#include "../src/io/dsl/file.h"
#include "../src/io/dsl/registry.h"
#include "../src/io/reader/commands/register_functions.h"
#include "../src/material/isotropic_j2_elasticity.h"
#include "../src/model/model.h"

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

using namespace fem;

TEST(Material_J2, DerivedModuliAndYieldPoints) {
    material::IsotropicJ2Elasticity j2(Precision(210000), Precision(0.3));

    EXPECT_NEAR(j2.shear_modulus(), Precision(80769.23076923077), Precision(1e-8));
    EXPECT_NEAR(j2.bulk_modulus(), Precision(175000), Precision(1e-8));

    j2.add_yield_point(Precision(280), Precision(0));
    j2.add_yield_point(Precision(320), Precision(0.02));
    j2.add_yield_point(Precision(380), Precision(0.10));

    const auto& points = j2.get_yield_points();
    ASSERT_EQ(points.size(), 3u);
    EXPECT_EQ(points[0].yield_stress, Precision(280));
    EXPECT_EQ(points[0].equivalent_plastic_strain, Precision(0));
    EXPECT_EQ(points[1].yield_stress, Precision(320));
    EXPECT_EQ(points[1].equivalent_plastic_strain, Precision(0.02));
    EXPECT_EQ(points[2].yield_stress, Precision(380));
    EXPECT_EQ(points[2].equivalent_plastic_strain, Precision(0.10));
}

TEST(Material_J2, PlasticCommandReplacesIsotropicElasticity) {
    const std::string input_path = "tests/TMP_J2_MATERIAL.INP";
    std::filesystem::remove(input_path);

    // Intentionally place PLASTIC before ELASTIC in the source. The semantic
    // material order must still execute ELASTIC first and PLASTIC afterwards.
    {
        std::ofstream os(input_path);
        ASSERT_TRUE(os.is_open());
        os << "*MATERIAL, NAME=STEEL\n";
        os << "*PLASTIC\n";
        os << "280., 0.\n";
        os << "320., 0.02\n";
        os << "380., 0.10\n";
        os << "*ELASTIC\n";
        os << "210000., 0.3\n";
    }

    model::Model model;
    io::dsl::Registry registry;
    io::reader::commands::register_material(registry, model);
    io::reader::commands::register_elastic(registry, model);
    io::reader::commands::register_plastic(registry, model);

    io::dsl::File file(input_path);
    io::dsl::DeckParser parser(registry);
    const io::dsl::Deck deck = parser.parse(file);

    const auto materials = deck.root().children("MATERIAL");
    ASSERT_EQ(materials.size(), 1u);

    const auto* definition = materials.front();
    definition->enter();
    definition->execute_children("ELASTIC");
    definition->execute_children("PLASTIC");
    definition->leave();

    const auto steel = model._data->materials.get("STEEL");
    ASSERT_NE(steel, nullptr);
    ASSERT_TRUE(steel->has_elasticity());

    const auto* j2 = steel->elasticity()->as<material::IsotropicJ2Elasticity>();
    ASSERT_NE(j2, nullptr);
    EXPECT_EQ(j2->youngs, Precision(210000));
    EXPECT_EQ(j2->poisson, Precision(0.3));

    const auto& points = j2->get_yield_points();
    ASSERT_EQ(points.size(), 3u);
    EXPECT_EQ(points[0].yield_stress, Precision(280));
    EXPECT_EQ(points[1].yield_stress, Precision(320));
    EXPECT_EQ(points[2].yield_stress, Precision(380));

    std::filesystem::remove(input_path);
}
