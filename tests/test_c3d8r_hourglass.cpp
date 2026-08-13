/**
 * @file test_c3d8r_hourglass.cpp
 * @brief Tests the parameter-free projected C3D8R hourglass stabilization.
 *
 * @author Finn Eggers
 * @date 13.08.2026
 */

#include "../src/material/isotropic_elasticity.h"
#include "../src/model/model.h"
#include "../src/model/solid/c3d8.h"
#include "../src/model/solid/c3d8r.h"

#include <gtest/gtest.h>

using namespace fem;

namespace {

template<class ElementType>
model::Model build_unit_cube(Precision youngs_modulus, Precision poisson_ratio) {
    model::Model model(8, 1, 0);

    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_node(2, 1.0, 1.0, 0.0);
    model.set_node(3, 0.0, 1.0, 0.0);
    model.set_node(4, 0.0, 0.0, 1.0);
    model.set_node(5, 1.0, 0.0, 1.0);
    model.set_node(6, 1.0, 1.0, 1.0);
    model.set_node(7, 0.0, 1.0, 1.0);
    model.set_element<ElementType>(0, 0, 1, 2, 3, 4, 5, 6, 7);

    auto material = model._data->materials.activate("MAT");
    material->set_elasticity<material::IsotropicElasticity>(youngs_modulus, poisson_ratio);
    model.solid_section("EALL", "MAT");
    model.assign_sections();

    return model;
}

StaticVector<24> affine_displacement() {
    static const Precision xyz[8][3] = {
        {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
        {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}
    };

    StaticVector<24> u = StaticVector<24>::Zero();
    for (Index node = 0; node < 8; ++node) {
        const Precision x = xyz[node][0];
        const Precision y = xyz[node][1];
        const Precision z = xyz[node][2];
        u(3 * node + 0) = Precision(0.10) * x + Precision(0.20) * y - Precision(0.05) * z;
        u(3 * node + 1) = Precision(0.03) * x - Precision(0.08) * y + Precision(0.04) * z;
        u(3 * node + 2) = Precision(0.02) * x + Precision(0.06) * y + Precision(0.12) * z;
    }
    return u;
}

StaticVector<24> shear_hourglass_displacement() {
    // u_x = s*t: zero center strain, trace-free fully integrated strain.
    static const Precision st[8] = {1, 1, -1, -1, -1, -1, 1, 1};

    StaticVector<24> u = StaticVector<24>::Zero();
    for (Index node = 0; node < 8; ++node) {
        u(3 * node) = st[node];
    }
    return u;
}

Precision quadratic_energy(const DynamicMatrix& stiffness, const StaticVector<24>& displacement) {
    return (displacement.transpose() * stiffness * displacement)(0, 0);
}

} // namespace

TEST(Elements_C3D8R, ProjectedHourglassDoesNotChangeAffineResponse) {
    auto full_model = build_unit_cube<model::C3D8>(1000.0, 0.3);
    auto reduced_model = build_unit_cube<model::C3D8R>(1000.0, 0.3);

    auto* full = full_model._data->elements[0]->as<model::C3D8>();
    auto* reduced = reduced_model._data->elements[0]->as<model::C3D8R>();
    ASSERT_NE(full, nullptr);
    ASSERT_NE(reduced, nullptr);

    Precision full_storage[24 * 24] {};
    Precision reduced_storage[24 * 24] {};
    const DynamicMatrix K_full = full->stiffness(full_storage);
    const DynamicMatrix K_reduced = reduced->stiffness(reduced_storage);
    const StaticVector<24> u = affine_displacement();

    const Precision full_energy = quadratic_energy(K_full, u);
    const Precision reduced_energy = quadratic_energy(K_reduced, u);
    EXPECT_NEAR(reduced_energy, full_energy,
                Precision(1e-10) * std::max(Precision(1), std::abs(full_energy)));
}

TEST(Elements_C3D8R, ProjectedHourglassRecoversPureShearHourglassEnergy) {
    auto full_model = build_unit_cube<model::C3D8>(1000.0, 0.3);
    auto reduced_model = build_unit_cube<model::C3D8R>(1000.0, 0.3);

    auto* full = full_model._data->elements[0]->as<model::C3D8>();
    auto* reduced = reduced_model._data->elements[0]->as<model::C3D8R>();
    ASSERT_NE(full, nullptr);
    ASSERT_NE(reduced, nullptr);

    Precision full_storage[24 * 24] {};
    Precision reduced_storage[24 * 24] {};
    const DynamicMatrix K_full = full->stiffness(full_storage);
    const DynamicMatrix K_reduced = reduced->stiffness(reduced_storage);
    const StaticVector<24> u = shear_hourglass_displacement();

    const Precision full_energy = quadratic_energy(K_full, u);
    const Precision reduced_energy = quadratic_energy(K_reduced, u);
    EXPECT_GT(full_energy, Precision(0));
    EXPECT_NEAR(reduced_energy, full_energy,
                Precision(1e-10) * std::max(Precision(1), std::abs(full_energy)));
}
