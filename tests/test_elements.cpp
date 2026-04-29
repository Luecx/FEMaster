#include "../src/material/abd_elasticity.h"
#include "../src/material/isotropic_elasticity.h"
#include "../src/model/model.h"
#include "../src/model/shell/qspt.h"
#include "../src/model/shell/s4.h"
#include "../src/model/truss/truss.h"

#include <gtest/gtest.h>

namespace {

fem::model::Model build_qspt_model(bool with_density) {
    fem::model::Model model(8, 4, 8);

    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_node(2, 1.0, 1.0, 0.0);
    model.set_node(3, 0.0, 1.0, 0.0);

    model.set_element<fem::model::QSPT>(0, 0, 1, 2, 3);

    auto material = model._data->materials.activate("MAT");
    material->set_elasticity<fem::material::IsotropicElasticity>(100.0, 0.0);
    if (with_density) {
        material->set_density(10.0);
    }

    model.shell_section("EALL", "MAT", 0.1);
    model.assign_sections();

    return model;
}

} // namespace

TEST(Elements_QSPT, StiffnessMassAndShearFlowForUnitSquare) {
    auto model = build_qspt_model(true);
    auto* elem = model._data->elements[0]->as<fem::model::QSPT>();
    ASSERT_NE(elem, nullptr);

    fem::Precision k_storage[12 * 12] {};
    fem::DynamicMatrix K = elem->stiffness(k_storage);
    EXPECT_TRUE(K.isApprox(K.transpose(), 1e-12));

    fem::StaticVector<12> u = fem::StaticVector<12>::Zero();
    u(6) = 1.0;
    u(9) = 1.0;
    const fem::Precision energy = (u.transpose() * K * u)(0, 0);
    EXPECT_NEAR(energy, 5.0, 1e-12);

    fem::Precision m_storage[12 * 12] {};
    fem::DynamicMatrix M = elem->mass(m_storage);
    EXPECT_TRUE(M.isApprox(M.transpose(), 1e-12));
    EXPECT_NEAR(M.sum(), 3.0, 1e-12);

    fem::model::Field displacement{"U", fem::model::FieldDomain::NODE, 4, 6};
    displacement.set_zero();
    displacement(2, 0) = 1.0;
    displacement(3, 0) = 1.0;

    auto shear_flow = model.compute_shear_flow(displacement);
    ASSERT_EQ(shear_flow.rows(), 4);
    ASSERT_EQ(shear_flow.cols(), 3);

    const fem::Precision expected[4] {-5.0, 5.0, -5.0, 5.0};
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(shear_flow(i, 0), 0.0, 1e-12);
        EXPECT_NEAR(shear_flow(i, 1), static_cast<fem::Precision>(i), 1e-12);
        EXPECT_NEAR(shear_flow(i, 2), expected[i], 1e-12);
    }

    auto compliance = model.compute_compliance(displacement);
    EXPECT_NEAR(compliance(0, 0), 5.0, 1e-12);
}

TEST(Elements_QSPT, MassMatrixIsZeroWithoutDensity) {
    auto model = build_qspt_model(false);
    auto* elem = model._data->elements[0]->as<fem::model::QSPT>();
    ASSERT_NE(elem, nullptr);

    fem::Precision storage[12 * 12] {};
    fem::DynamicMatrix M = elem->mass(storage);
    EXPECT_NEAR(M.norm(), 0.0, 1e-12);
}

TEST(Elements_S4, ABDMaterialUsesMaterialDensityForMass) {
    fem::model::Model model(8, 4, 8);

    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_node(2, 1.0, 1.0, 0.0);
    model.set_node(3, 0.0, 1.0, 0.0);
    model.set_element<fem::model::S4>(0, 0, 1, 2, 3);

    auto material = model._data->materials.activate("MAT");
    material->set_density(10.0);

    fem::StaticMatrix<6, 6> abd = fem::StaticMatrix<6, 6>::Identity();
    fem::StaticMatrix<2, 2> shear = fem::StaticMatrix<2, 2>::Identity();
    material->set_elasticity<fem::material::ABDElasticity>(abd, shear);

    model.shell_section("EALL", "MAT", 0.1);
    model.assign_sections();

    auto* elem = model._data->elements[0]->as<fem::model::S4>();
    ASSERT_NE(elem, nullptr);

    fem::Precision k_storage[24 * 24] {};
    fem::DynamicMatrix K = elem->stiffness(k_storage);
    EXPECT_TRUE(K.isApprox(K.transpose(), 1e-12));

    fem::Precision m_storage[24 * 24] {};
    fem::DynamicMatrix M = elem->mass(m_storage);
    EXPECT_TRUE(M.isApprox(M.transpose(), 1e-12));

    fem::Precision ux_mass = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            ux_mass += M(6 * i, 6 * j);
        }
    }
    EXPECT_NEAR(ux_mass, 1.0, 1e-12);
}

TEST(Elements_Truss, UsesDedicatedTrussSectionArea) {
    fem::model::Model model(4, 4, 4);

    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_element<fem::model::T3>(0, 0, 1);

    auto material = model._data->materials.activate("MAT");
    material->set_elasticity<fem::material::IsotropicElasticity>(6.0, 0.0);
    material->set_density(3.0);

    model.truss_section("EALL", "MAT", 2.0);
    model.assign_sections();

    auto* elem = model._data->elements[0]->as<fem::model::T3>();
    ASSERT_NE(elem, nullptr);
    ASSERT_NE(elem->get_section(), nullptr);
    EXPECT_NEAR(elem->get_section()->A, 2.0, 1e-12);

    fem::Precision k_storage[6 * 6] {};
    fem::DynamicMatrix K = elem->stiffness(k_storage);
    EXPECT_NEAR(K(0, 0), 12.0, 1e-12);
    EXPECT_NEAR(K(0, 3), -12.0, 1e-12);
    EXPECT_NEAR(K(3, 0), -12.0, 1e-12);
    EXPECT_NEAR(K(3, 3), 12.0, 1e-12);

    fem::Precision m_storage[6 * 6] {};
    fem::DynamicMatrix M = elem->mass(m_storage);
    EXPECT_NEAR(M(0, 0), 3.0, 1e-12);
    EXPECT_NEAR(M(3, 3), 3.0, 1e-12);

    fem::model::Field displacement{"U", fem::model::FieldDomain::NODE, 2, 6};
    displacement.set_zero();
    displacement(1, 0) = 1.0;

    const auto section_forces = model.compute_section_forces(displacement);
    ASSERT_EQ(section_forces.rows(), 2);
    ASSERT_EQ(section_forces.cols(), 8);
    EXPECT_NEAR(section_forces(0, 0), 0.0, 1e-12);
    EXPECT_NEAR(section_forces(1, 0), 0.0, 1e-12);
    EXPECT_NEAR(section_forces(0, 1), 0.0, 1e-12);
    EXPECT_NEAR(section_forces(1, 1), 1.0, 1e-12);
    EXPECT_NEAR(section_forces(0, 2), 12.0, 1e-12);
    EXPECT_NEAR(section_forces(1, 2), 12.0, 1e-12);
    for (int col = 3; col < 8; ++col) {
        EXPECT_NEAR(section_forces(0, col), 0.0, 1e-12);
        EXPECT_NEAR(section_forces(1, col), 0.0, 1e-12);
    }
}

