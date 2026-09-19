/**
 * @file test_model_thermal.cpp
 * @brief Tests model-level thermal system helpers.
 */

#include "../src/material/isotropic_elasticity.h"
#include "../src/model/model.h"
#include "../src/model/solid/c3d8.h"
#include "../src/model/solid/c3d8r.h"
#include "../src/model/truss/truss.h"
#include "../src/section/section_solid.h"

#include <gtest/gtest.h>

using namespace fem;

TEST(ModelThermalDofs, ActivatesOnlyThermalElementNodes) {
    model::Model model;

    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_node(2, 1.0, 1.0, 0.0);
    model.set_node(3, 0.0, 1.0, 0.0);
    model.set_node(4, 0.0, 0.0, 1.0);
    model.set_node(5, 1.0, 0.0, 1.0);
    model.set_node(6, 1.0, 1.0, 1.0);
    model.set_node(7, 0.0, 1.0, 1.0);
    model.set_node(8, 2.0, 0.0, 0.0);
    model.set_node(9, 3.0, 0.0, 0.0);

    model.set_element<model::C3D8>(0, 0, 1, 2, 3, 4, 5, 6, 7);
    model.set_element<model::T3>(1, 8, 9);
    model.compile();

    const SystemDofIds ids = model.build_thermal_dof_index_matrix();

    ASSERT_EQ(ids.rows(), 10);
    ASSERT_EQ(ids.cols(), 1);

    for (Index node = 0; node < 8; ++node) {
        EXPECT_EQ(ids(node, 0), node);
    }

    EXPECT_EQ(ids(8, 0), -1);
    EXPECT_EQ(ids(9, 0), -1);
}


TEST(ModelThermalResults, RecoversReducedSolidHeatFluxAsNodalField) {
    model::Model model;

    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_node(2, 1.0, 1.0, 0.0);
    model.set_node(3, 0.0, 1.0, 0.0);
    model.set_node(4, 0.0, 0.0, 1.0);
    model.set_node(5, 1.0, 0.0, 1.0);
    model.set_node(6, 1.0, 1.0, 1.0);
    model.set_node(7, 0.0, 1.0, 1.0);
    model.set_element<model::C3D8R>(0, 0, 1, 2, 3, 4, 5, 6, 7);

    // A linear temperature field T=x has the exact constant gradient [1,0,0].
    // With k=2, Fourier's law therefore gives q=[-2,0,0] everywhere.
    auto material = std::make_shared<material::Material>("THERMAL");
    material->set_thermal_conductivity(Precision(2));
    model.add_material(material);

    const auto part = model._data->parts.get();
    ASSERT_NE(part, nullptr);

    auto section = std::make_shared<SolidSection>();
    section->material_ = material;
    section->region_   = part->elem_sets.get(SET_ELEM_ALL);
    model.add_section(section);

    model.compile();
    model.assign_sections();

    model::Field temperature{"TEMPERATURE", model::FieldDomain::NODE, 8, 1};
    for (Index node = 0; node < 8; ++node) {
        temperature(node, 0) = (node == 1 || node == 2 || node == 5 || node == 6)
            ? Precision(1)
            : Precision(0);
    }

    // C3D8R evaluates heat flux at its single integration point. The solid
    // formulation must extrapolate that value to ELEMENT_NODAL storage before
    // the model projects it to unique global nodes.
    const model::Field heat_flux = model.compute_heat_flux(temperature);

    EXPECT_EQ(heat_flux.domain, model::FieldDomain::NODE);
    ASSERT_EQ(heat_flux.rows, 8);
    ASSERT_EQ(heat_flux.components, 3);

    for (Index node = 0; node < heat_flux.rows; ++node) {
        EXPECT_NEAR(heat_flux(node, 0), Precision(-2), Precision(1e-12));
        EXPECT_NEAR(heat_flux(node, 1), Precision(0), Precision(1e-12));
        EXPECT_NEAR(heat_flux(node, 2), Precision(0), Precision(1e-12));
    }
}
