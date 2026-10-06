/**
 * @file test_solid_thermal_recovery.cpp
 * @brief Thermal free strain recovery and prestress for linear solid mechanics.
 */

#include "../src/material/isotropic_elasticity.h"
#include "../src/model/model.h"
#include "../src/model/solid/c3d8.h"
#include "../src/section/section_solid.h"

#include <gtest/gtest.h>

#include <array>
#include <memory>
#include <tuple>

using namespace fem;

TEST(SolidThermal, FreeExpansionAndRestrainedPrestress) {
    model::Model model;

    const std::array<Vec3, 8> coords{
        Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(1, 1, 0), Vec3(0, 1, 0),
        Vec3(0, 0, 1), Vec3(1, 0, 1), Vec3(1, 1, 1), Vec3(0, 1, 1)
    };

    for (Index node = 0; node < 8; ++node) {
        model.set_node(node, coords[node].x(), coords[node].y(), coords[node].z());
    }
    model.set_element<model::C3D8>(0, 0, 1, 2, 3, 4, 5, 6, 7);

    auto material = std::make_shared<material::Material>("MAT");
    material->set_elasticity<material::IsotropicElasticity>(1000.0, 0.25);
    material->set_thermal_expansion(0.01);
    material->set_thermal_zero_temperature(20.0);
    model.add_material(material);

    auto section = std::make_shared<SolidSection>();
    section->material_ = material;
    section->region_ = model._data->parts.get()->elem_sets.get(SET_ELEM_ALL);
    model.add_section(section);

    model.compile();
    model.step_begin();

    auto temperature = std::make_shared<model::Field>(
        "TEMP", model::FieldDomain::NODE, 8, 1);
    for (Index node = 0; node < 8; ++node) {
        (*temperature)(node, 0) = 40.0;
    }

    model._data->temperature = temperature;

    model::Field displacement{"DISPLACEMENT", model::FieldDomain::NODE, 8, 6};
    displacement.set_zero();
    for (Index node = 0; node < 8; ++node) {
        displacement(node, 0) = 0.2 * coords[node].x();
        displacement(node, 1) = 0.2 * coords[node].y();
        displacement(node, 2) = 0.2 * coords[node].z();
    }

    auto stress_strain =
        model.compute_stress_nodal(displacement);
    const auto& stress = std::get<0>(stress_strain);
    const auto& strain = std::get<1>(stress_strain);

    for (Index node = 0; node < 8; ++node) {
        EXPECT_NEAR(strain(node, 0), 0.2, 1e-9);
        EXPECT_NEAR(strain(node, 1), 0.2, 1e-9);
        EXPECT_NEAR(strain(node, 2), 0.2, 1e-9);
        for (Index component = 0; component < 6; ++component) {
            EXPECT_NEAR(stress(node, component), 0.0, 1e-8)
                << "node=" << node << ", component=" << component;
        }
    }

    auto* element = model._data->elements[0]->as<model::C3D8>();
    ASSERT_NE(element, nullptr);

    Precision kg_storage[24 * 24] {};
    const DynamicMatrix Kg_free =
        element->evaluate(nullptr, kg_storage, nullptr, &displacement, nullptr, false);
    EXPECT_LT(Kg_free.norm(), 1e-8);

    displacement.set_zero();
    stress_strain =
        model.compute_stress_nodal(displacement);
    const auto& restrained_stress = std::get<0>(stress_strain);
    const Precision expected = -1000.0 / (1.0 - 2.0 * 0.25) * 0.2;

    for (Index node = 0; node < 8; ++node) {
        EXPECT_NEAR(restrained_stress(node, 0), expected, 1e-8);
        EXPECT_NEAR(restrained_stress(node, 1), expected, 1e-8);
        EXPECT_NEAR(restrained_stress(node, 2), expected, 1e-8);
        EXPECT_NEAR(restrained_stress(node, 3), 0.0, 1e-8);
        EXPECT_NEAR(restrained_stress(node, 4), 0.0, 1e-8);
        EXPECT_NEAR(restrained_stress(node, 5), 0.0, 1e-8);
    }

    const DynamicMatrix Kg_restrained =
        element->evaluate(nullptr, kg_storage, nullptr, &displacement, nullptr, false);
    EXPECT_GT(Kg_restrained.norm(), 1e-6);

    model.step_end();
}
