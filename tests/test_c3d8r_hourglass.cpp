/**
 * @file test_c3d8r_hourglass.cpp
 * @brief Tests the parameter-free physical C3D8R hourglass stabilization.
 *
 * @author Finn Eggers
 * @date 13.08.2026
 */

#include "../src/material/isotropic_elasticity.h"
#include "../src/model/model.h"
#include "../src/model/solid/c3d8.h"
#include "../src/model/solid/c3d8r.h"

#include <algorithm>
#include <cmath>
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

model::Model build_rotated_unit_cube(Precision youngs_modulus, Precision poisson_ratio) {
    model::Model model(8, 1, 0);
    static const Precision xyz[8][3] = {
        {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
        {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}
    };

    // Rigid +90 degree rotation about global z followed by translation.
    for (Index node = 0; node < 8; ++node) {
        const Precision x = xyz[node][0];
        const Precision y = xyz[node][1];
        const Precision z = xyz[node][2];
        model.set_node(node, Precision(2) - y, Precision(3) + x, Precision(4) + z);
    }

    model.set_element<model::C3D8R>(0, 0, 1, 2, 3, 4, 5, 6, 7);
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

StaticVector<24> shear_hourglass_displacement(Dim direction = 0) {
    static const Precision st[8] = {1, 1, -1, -1, -1, -1, 1, 1};

    StaticVector<24> u = StaticVector<24>::Zero();
    for (Index node = 0; node < 8; ++node) {
        u(3 * node + direction) = st[node];
    }
    return u;
}

Precision quadratic_energy(const DynamicMatrix& stiffness, const StaticVector<24>& displacement) {
    return (displacement.transpose() * stiffness * displacement)(0, 0);
}

} // namespace

TEST(Elements_C3D8R, PhysicalHourglassDoesNotChangeAffineResponse) {
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

TEST(Elements_C3D8R, PhysicalHourglassMatchesUnitCubeAssumedStrainEnergy) {
    constexpr Precision youngs_modulus = 1000.0;
    constexpr Precision poisson_ratio  = 0.3;
    const Precision shear_modulus =
        youngs_modulus / (Precision(2) * (Precision(1) + poisson_ratio));

    auto reduced_model = build_unit_cube<model::C3D8R>(youngs_modulus, poisson_ratio);
    auto* reduced = reduced_model._data->elements[0]->as<model::C3D8R>();
    ASSERT_NE(reduced, nullptr);

    Precision reduced_storage[24 * 24] {};
    const DynamicMatrix K_reduced = reduced->stiffness(reduced_storage);
    const StaticVector<24> u = shear_hourglass_displacement();

    // Unit cube: H_11 = H_22 = H_33 = 4/3 and gamma_1^T gamma_1 = 8.
    // For u_x = gamma_1, u^T K_hg u = 2 mu (4/3) (8)^2.
    const Precision expected_energy =
        Precision(512) / Precision(3) * shear_modulus;
    const Precision reduced_energy = quadratic_energy(K_reduced, u);

    EXPECT_GT(reduced_energy, Precision(0));
    EXPECT_NEAR(reduced_energy, expected_energy,
                Precision(1e-10) * std::max(Precision(1), std::abs(expected_energy)));
}

TEST(Elements_C3D8R, PhysicalHourglassIsInvariantToRigidElementOrientation) {
    constexpr Precision youngs_modulus = 1000.0;
    constexpr Precision poisson_ratio  = 0.3;

    auto reference_model = build_unit_cube<model::C3D8R>(youngs_modulus, poisson_ratio);
    auto rotated_model   = build_rotated_unit_cube(youngs_modulus, poisson_ratio);
    auto* reference = reference_model._data->elements[0]->as<model::C3D8R>();
    auto* rotated   = rotated_model._data->elements[0]->as<model::C3D8R>();
    ASSERT_NE(reference, nullptr);
    ASSERT_NE(rotated, nullptr);

    Precision reference_storage[24 * 24] {};
    Precision rotated_storage[24 * 24] {};
    const DynamicMatrix K_reference = reference->stiffness(reference_storage);
    const DynamicMatrix K_rotated   = rotated->stiffness(rotated_storage);

    // Local x rotates to global y under the +90 degree rigid rotation.
    const Precision reference_energy =
        quadratic_energy(K_reference, shear_hourglass_displacement(0));
    const Precision rotated_energy =
        quadratic_energy(K_rotated, shear_hourglass_displacement(1));

    EXPECT_NEAR(rotated_energy, reference_energy,
                Precision(1e-10) * std::max(Precision(1), std::abs(reference_energy)));
}

TEST(Elements_C3D8R, PhysicalHourglassRemainsBoundedNearIncompressibility) {
    constexpr Precision shear_modulus = 1000.0;
    constexpr Precision poisson_ratio = 0.4999;
    const Precision youngs_modulus =
        Precision(2) * shear_modulus * (Precision(1) + poisson_ratio);

    auto reduced_model = build_unit_cube<model::C3D8R>(youngs_modulus, poisson_ratio);
    auto* reduced = reduced_model._data->elements[0]->as<model::C3D8R>();
    ASSERT_NE(reduced, nullptr);

    Precision reduced_storage[24 * 24] {};
    const DynamicMatrix K_reduced = reduced->stiffness(reduced_storage);
    const StaticVector<24> u = shear_hourglass_displacement();

    const Precision expected_energy =
        Precision(512) / Precision(3) * shear_modulus;
    const Precision reduced_energy = quadratic_energy(K_reduced, u);

    EXPECT_TRUE(std::isfinite(reduced_energy));
    EXPECT_NEAR(reduced_energy, expected_energy,
                Precision(1e-9) * std::max(Precision(1), std::abs(expected_energy)));
}
