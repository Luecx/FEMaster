/**
 * @file test_c3d8i.cpp
 * @brief Regression tests for the incompatible-mode C3D8I solid.
 */

#include "../src/loadcase/tools/nonlinear_state_manager.h"
#include "../src/material/isotropic_elasticity.h"
#include "../src/material/neo_hooke_elasticity.h"
#include "../src/model/model.h"
#include "../src/model/solid/c3d8.h"
#include "../src/model/solid/c3d8i.h"
#include "../src/section/section_solid.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <memory>

using namespace fem;

namespace {

template<class Element, class Elasticity, class... Args>
void build_unit_hex(model::Model& model, Args&&... args) {
    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_node(2, 1.0, 1.0, 0.0);
    model.set_node(3, 0.0, 1.0, 0.0);
    model.set_node(4, 0.0, 0.0, 1.0);
    model.set_node(5, 1.0, 0.0, 1.0);
    model.set_node(6, 1.0, 1.0, 1.0);
    model.set_node(7, 0.0, 1.0, 1.0);
    model.set_element<Element>(0, 0, 1, 2, 3, 4, 5, 6, 7);

    auto material = std::make_shared<material::Material>("MAT");
    material->set_elasticity<Elasticity>(std::forward<Args>(args)...);
    model.add_material(material);

    auto section = std::make_shared<SolidSection>();
    section->material_ = material;
    section->region_ = model._data->parts.get()->elem_sets.get(SET_ELEM_ALL);
    model.add_section(section);

    model.compile();
    model.assign_sections();
    model.step_begin();
}

StaticVector<24> local_force(const model::NodeData& forces) {
    StaticVector<24> result = StaticVector<24>::Zero();
    for (Index a = 0; a < 8; ++a) {
        for (Dim d = 0; d < 3; ++d) result(3 * a + d) = forces(a, d);
    }
    return result;
}

model::Field zero_displacement() {
    model::Field displacement("U", model::FieldDomain::NODE, 8, 6);
    displacement.set_zero();
    return displacement;
}

} // namespace

TEST(Elements_C3D8I, AffineLinearPatchMatchesC3D8) {
    model::Model compatible;
    model::Model incompatible;

    build_unit_hex<model::C3D8, material::IsotropicElasticity>(
        compatible, 1000.0, 0.3);
    build_unit_hex<model::C3D8I, material::IsotropicElasticity>(
        incompatible, 1000.0, 0.3);

    std::array<Precision, 24 * 24> storage_c {};
    std::array<Precision, 24 * 24> storage_i {};

    const Matrix24 Kc = compatible._data->elements[0]
        ->as<model::C3D8>()->stiffness(storage_c.data());
    const Matrix24 Ki = incompatible._data->elements[0]
        ->as<model::C3D8I>()->stiffness(storage_i.data());

    StaticVector<24> u = StaticVector<24>::Zero();
    for (Index a = 0; a < 8; ++a) {
        const Vec3 X = incompatible._data->positions_reference->row_vec3(a);
        u(3 * a + 0) = 0.013 * X(0) + 0.004 * X(1) - 0.002 * X(2);
        u(3 * a + 1) = 0.003 * X(0) - 0.007 * X(1) + 0.005 * X(2);
        u(3 * a + 2) = 0.002 * X(0) + 0.006 * X(1) + 0.009 * X(2);
    }

    const auto fc = Kc * u;
    const auto fi = Ki * u;

    EXPECT_LT((fc - fi).norm() / std::max(Precision(1), fc.norm()), 1e-11);

    compatible.step_end();
    incompatible.step_end();
}

TEST(Elements_C3D8I, FiniteRigidRotationIsObjective) {
    model::Model model;
    build_unit_hex<model::C3D8I, material::NeoHookeElasticity>(
        model, 0.5, 0.5);

    loadcase::tools::NonlinearStateManager state(model);
    state.begin_element_trial();

    auto displacement = zero_displacement();
    const Precision angle = 0.63;
    Mat3 R = Mat3::Identity();
    R(0, 0) =  std::cos(angle);
    R(0, 1) = -std::sin(angle);
    R(1, 0) =  std::sin(angle);
    R(1, 1) =  std::cos(angle);

    for (Index a = 0; a < 8; ++a) {
        const Vec3 X = model._data->positions_reference->row_vec3(a);
        const Vec3 u = R * X - X;
        displacement(a, 0) = u(0);
        displacement(a, 1) = u(1);
        displacement(a, 2) = u(2);
    }

    model::NodeData forces(
        "F", model::FieldDomain::NODE, 8, 6);
    forces.set_zero();

    std::array<Precision, 24 * 24> storage {};
    model._data->elements[0]->as<model::C3D8I>()
        ->stiffness_tangent(storage.data(), forces, displacement);

    EXPECT_LT(local_force(forces).norm(), 1e-9);

    state.rollback_element_trial();
    model.step_end();
}

TEST(Elements_C3D8I, NonlinearCondensedTangentMatchesFiniteDifference) {
    model::Model model;
    build_unit_hex<model::C3D8I, material::NeoHookeElasticity>(
        model, 0.7, 0.8);

    loadcase::tools::NonlinearStateManager state(model);
    state.begin_element_trial();

    auto displacement = zero_displacement();
    for (Index a = 0; a < 8; ++a) {
        const Vec3 X = model._data->positions_reference->row_vec3(a);
        displacement(a, 0) =
            0.08 * X(0) + 0.035 * X(0) * X(1) - 0.015 * X(2);
        displacement(a, 1) =
           -0.03 * X(1) + 0.025 * X(1) * X(2) + 0.010 * X(0);
        displacement(a, 2) =
            0.05 * X(2) + 0.020 * X(0) * X(2) - 0.012 * X(1);
    }

    auto* element = model._data->elements[0]->as<model::C3D8I>();

    model::NodeData forces(
        "F", model::FieldDomain::NODE, 8, 6);
    forces.set_zero();

    std::array<Precision, 24 * 24> storage {};
    const Matrix24 tangent =
        element->stiffness_tangent(storage.data(), forces, displacement);

    Matrix24 finite_difference = Matrix24::Zero();
    const Precision h = 2e-7;

    for (Index column = 0; column < 24; ++column) {
        const Index node = column / 3;
        const Index dof  = column % 3;

        auto plus = displacement;
        auto minus = displacement;
        plus(node, dof)  += h;
        minus(node, dof) -= h;

        model::NodeData f_plus(
            "F+", model::FieldDomain::NODE, 8, 6);
        model::NodeData f_minus(
            "F-", model::FieldDomain::NODE, 8, 6);
        f_plus.set_zero();
        f_minus.set_zero();

        element->stiffness_tangent(nullptr, f_plus, plus);
        element->stiffness_tangent(nullptr, f_minus, minus);

        finite_difference.col(column) =
            (local_force(f_plus) - local_force(f_minus)) / (2 * h);
    }

    const Precision relative_error =
        (tangent - finite_difference).norm()
        / std::max(Precision(1), finite_difference.norm());

    EXPECT_LT(relative_error, 2e-5);

    state.rollback_element_trial();
    model.step_end();
}
