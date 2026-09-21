/**
 * @file test_shell_thermal.cpp
 * @brief Uniform-through-thickness thermal equivalent loads for MITC/FRT shells.
 */

#include "../src/bc/neumann/load_t.h"
#include "../src/material/isotropic_elasticity.h"
#include "../src/model/model.h"
#include "../src/model/shell/frt_shell_s3.h"
#include "../src/model/shell/frt_shell_s4.h"
#include "../src/model/shell/frt_shell_s6.h"
#include "../src/model/shell/frt_shell_s8.h"
#include "../src/section/section_shell_abd.h"
#include "../src/section/section_shell_integrated.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <memory>
#include <utility>

namespace {

template<typename Shell, std::size_t... I>
void add_shell(fem::model::Model& model, std::index_sequence<I...>) {
    model.set_element<Shell>(0, static_cast<fem::ID>(I)...);
}

template<typename Shell, std::size_t N>
void check_uniform_thermal_load(
    const std::array<fem::Vec3, N>& coords, bool use_abd
) {
    using namespace fem;

    model::Model model;
    for (std::size_t i = 0; i < N; ++i) {
        model.set_node(static_cast<ID>(i), coords[i].x(), coords[i].y(), coords[i].z());
    }
    add_shell<Shell>(model, std::make_index_sequence<N>{});

    auto material = std::make_shared<material::Material>("MAT");
    material->set_elasticity<material::IsotropicElasticity>(1000.0, 0.25);
    material->set_thermal_expansion(0.01);
    model.add_material(material);

    const auto region = model._data->parts.get()->elem_sets.get(SET_ELEM_ALL);
    if (use_abd) {
        // Symmetric positive-definite ABD with explicit membrane/bending coupling.
        Mat6 abd = Mat6::Zero();
        abd.diagonal() << 100.0, 100.0, 40.0, 10.0, 10.0, 4.0;
        abd(0, 3) = abd(3, 0) = 5.0;
        abd(1, 4) = abd(4, 1) = 5.0;
        Mat2 shear = 20.0 * Mat2::Identity();
        model.add_section(std::make_shared<ABDShellSection>(
            material, region, 0.2, abd, shear, nullptr
        ));
    } else {
        model.add_section(std::make_shared<IntegratedShellSection>(
            material, region, 0.2, nullptr
        ));
    }

    model.compile();
    model.step_begin();

    auto* element = model._data->elements[0]->as<Shell>();
    ASSERT_NE(element, nullptr);

    auto temperature = std::make_shared<model::Field>(
        "TEMP", model::FieldDomain::NODE, static_cast<Index>(N), 1
    );
    for (Index i = 0; i < static_cast<Index>(N); ++i) {
        (*temperature)(i, 0) = 40.0;
    }

    bc::TLoad load;
    load.temp_field_ = temperature;
    load.ref_temp_ = 20.0;

    model::Field rhs{"RHS", model::FieldDomain::NODE, static_cast<Index>(N), 6};
    rhs.set_zero();
    load.apply(*model._data, rhs, 0.0);

    // A homogeneous thermal free-expansion displacement produces exactly the
    // same generalized strain as the uniform, thickness-constant thermal load.
    // Comparing with K*u checks all six DOFs, correct sign, thickness factor,
    // topology-specific MITC B, and (for ABD) coupled nodal thermal moments.
    Precision matrix_storage[6 * N * 6 * N] {};
    const DynamicMatrix K = element->stiffness(matrix_storage);
    StaticVector<6 * N> free_expansion = StaticVector<6 * N>::Zero();
    for (Index node = 0; node < static_cast<Index>(N); ++node) {
        free_expansion(6 * node + 0) = 0.2 * coords[node].x();
        free_expansion(6 * node + 1) = 0.2 * coords[node].y();
    }
    const auto expected = (K * free_expansion).eval();

    Precision rotational_load = 0.0;
    for (Index node = 0; node < static_cast<Index>(N); ++node) {
        for (Index dof = 0; dof < 6; ++dof) {
            EXPECT_NEAR(rhs(node, dof), expected(6 * node + dof), 1e-8)
                << "node=" << node << ", dof=" << dof;
        }
        rotational_load += std::abs(rhs(node, 3)) + std::abs(rhs(node, 4));
    }
    if (use_abd) {
        EXPECT_GT(rotational_load, 1e-6);
    }

    // No thermal RHS when all nodal temperatures equal the reference state.
    for (Index node = 0; node < static_cast<Index>(N); ++node) {
        (*temperature)(node, 0) = load.ref_temp_;
    }
    rhs.set_zero();
    load.apply(*model._data, rhs, 0.0);
    for (Index node = 0; node < static_cast<Index>(N); ++node) {
        for (Index dof = 0; dof < 6; ++dof) {
            EXPECT_NEAR(rhs(node, dof), 0.0, 1e-12);
        }
    }

    model.step_end();
}

} // namespace

TEST(ShellThermal, S3Integrated) {
    check_uniform_thermal_load<fem::model::FRTShellS3>(
        {fem::Vec3(0, 0, 0), fem::Vec3(1, 0, 0), fem::Vec3(0, 1, 0)}, false
    );
}

TEST(ShellThermal, S4Integrated) {
    check_uniform_thermal_load<fem::model::FRTShellS4>(
        {fem::Vec3(0, 0, 0), fem::Vec3(1, 0, 0),
         fem::Vec3(1, 1, 0), fem::Vec3(0, 1, 0)}, false
    );
}

TEST(ShellThermal, S6Integrated) {
    check_uniform_thermal_load<fem::model::FRTShellS6>(
        {fem::Vec3(0, 0, 0), fem::Vec3(1, 0, 0), fem::Vec3(0, 1, 0),
         fem::Vec3(0.5, 0, 0), fem::Vec3(0.5, 0.5, 0),
         fem::Vec3(0, 0.5, 0)}, false
    );
}

TEST(ShellThermal, S8Integrated) {
    check_uniform_thermal_load<fem::model::FRTShellS8>(
        {fem::Vec3(0, 0, 0), fem::Vec3(1, 0, 0),
         fem::Vec3(1, 1, 0), fem::Vec3(0, 1, 0),
         fem::Vec3(0.5, 0, 0), fem::Vec3(1, 0.5, 0),
         fem::Vec3(0.5, 1, 0), fem::Vec3(0, 0.5, 0)}, false
    );
}

TEST(ShellThermal, S4ABDWithMembraneBendingCoupling) {
    check_uniform_thermal_load<fem::model::FRTShellS4>(
        {fem::Vec3(0, 0, 0), fem::Vec3(1, 0, 0),
         fem::Vec3(1, 1, 0), fem::Vec3(0, 1, 0)}, true
    );
}
