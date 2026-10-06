/**
 * @file test_b31.cpp
 * @brief Tests the geometrically exact two-node B31 beam.
 */

#include "../src/material/isotropic_elasticity.h"
#include "../src/math/so3.h"
#include "../src/model/beam/b31.h"
#include "../src/model/model.h"
#include "../src/section/profile.h"
#include "../src/section/section_beam.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <memory>

using namespace fem;

namespace {

model::Model build_b31_model(bool with_offsets = true) {
    model::Model model;

    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 2.0, 0.0, 0.0);
    model.set_element<model::B31>(0, 0, 1);

    auto material = std::make_shared<material::Material>("MAT");
    material->set_elasticity<material::IsotropicElasticity>(1000.0, 0.25);
    material->set_density(3.0);
    model.add_material(material);

    auto profile = std::make_shared<Profile>(
        "P",
        2.0,
        0.30,
        0.40,
        0.10,
        0.02,
        with_offsets ? 0.12 : 0.0,
        with_offsets ? -0.08 : 0.0,
        with_offsets ? -0.05 : 0.0,
        with_offsets ? 0.06 : 0.0,
        1.50,
        1.40
    );
    model.add_profile(profile);

    auto section = std::make_shared<BeamSection>();
    section->material_  = material;
    section->profile_   = profile;
    section->direction_ = Vec3(0.0, 1.0, 0.0);
    section->region_    = model._data->parts.get()->elem_sets.get(SET_ELEM_ALL);
    model.add_section(section);

    model.compile();
    return model;
}

model::NodeData exact_force(model::B31& element, const model::Field& displacement) {
    model::NodeData force("F", model::FieldDomain::NODE, 2, 6);
    force.set_zero();
    element.evaluate(
        nullptr,
        nullptr,
        &force,
        &displacement,
        nullptr,
        &displacement,
        nullptr,
        false
    );
    return force;
}

StaticVector<12> gather_force(const model::NodeData& force) {
    StaticVector<12> result = StaticVector<12>::Zero();
    for (Index node = 0; node < 2; ++node) {
        for (Index dof = 0; dof < 6; ++dof) {
            result(6 * node + dof) = force(node, dof);
        }
    }
    return result;
}

} // namespace

TEST(Elements_B31, ReferenceTangentAndMassAreSymmetric) {
    auto model = build_b31_model();
    auto* element = model._data->elements[0]->as<model::B31>();
    ASSERT_NE(element, nullptr);

    Precision stiffness_storage[12 * 12] {};
    DynamicMatrix K = element->evaluate(
        stiffness_storage,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        false
    );

    Precision mass_storage[12 * 12] {};
    DynamicMatrix M = element->mass(mass_storage);

    EXPECT_TRUE(K.allFinite());
    EXPECT_TRUE(M.allFinite());
    EXPECT_TRUE(K.isApprox(K.transpose(), 1e-10));
    EXPECT_TRUE(M.isApprox(M.transpose(), 1e-12));
    EXPECT_GT(K.norm(), Precision(0));
    EXPECT_GT(M.norm(), Precision(0));
}

TEST(Elements_B31, FiniteRigidBodyMotionIsExactlyStrainFreeWithOffsets) {
    auto model = build_b31_model(true);
    auto* element = model._data->elements[0]->as<model::B31>();
    ASSERT_NE(element, nullptr);

    const Vec3 theta(0.45, -0.25, 0.30);
    const Mat3 R = math::so3::rotation_matrix(theta);
    const Vec3 translation(0.35, -0.20, 0.15);

    model::Field displacement("U", model::FieldDomain::NODE, 2, 6);
    displacement.set_zero();

    const std::array<Vec3, 2> X{
        Vec3(0.0, 0.0, 0.0),
        Vec3(2.0, 0.0, 0.0)
    };

    for (Index node = 0; node < 2; ++node) {
        const Vec3 u = translation + R * X[node] - X[node];
        for (Index d = 0; d < 3; ++d) {
            displacement(node, d)     = u(d);
            displacement(node, d + 3) = theta(d);
        }
    }

    const auto force = exact_force(*element, displacement);
    EXPECT_LT(gather_force(force).norm(), Precision(1e-8));
}

TEST(Elements_B31, NonlinearTangentMatchesExactResidualDerivative) {
    auto model = build_b31_model(true);
    auto* element = model._data->elements[0]->as<model::B31>();
    ASSERT_NE(element, nullptr);

    model::Field state("U", model::FieldDomain::NODE, 2, 6);
    state.set_zero();
    state(0, 1) = -0.03;
    state(0, 3) =  0.08;
    state(1, 0) =  0.06;
    state(1, 1) =  0.15;
    state(1, 2) = -0.04;
    state(1, 3) =  0.12;
    state(1, 4) = -0.18;
    state(1, 5) =  0.22;

    Precision storage[12 * 12] {};
    DynamicMatrix K = element->evaluate(
        storage,
        nullptr,
        nullptr,
        &state,
        nullptr,
        &state,
        nullptr,
        false
    );

    const Precision h_translation = Precision(2e-6);
    const Precision h_rotation    = Precision(2e-6);

    for (Index column = 0; column < 12; ++column) {
        model::Field plus = state;
        model::Field minus = state;

        const Precision h = (column % 6) < 3
            ? h_translation
            : h_rotation;

        plus(column / 6, column % 6)  += h;
        minus(column / 6, column % 6) -= h;

        const StaticVector<12> f_plus =
            gather_force(exact_force(*element, plus));
        const StaticVector<12> f_minus =
            gather_force(exact_force(*element, minus));
        const StaticVector<12> numerical =
            (f_plus - f_minus) / (Precision(2) * h);

        const Precision scale =
            std::max(Precision(1), numerical.norm());
        EXPECT_LT(
            (K.col(column) - numerical).norm() / scale,
            Precision(2e-4)
        ) << "column " << column;
    }
}

TEST(Elements_B31, DefaultShearAreasAreFiveSixthsOfArea) {
    Profile profile("P", 12.0, 2.0, 3.0, 1.0);
    EXPECT_NEAR(profile.shear_area_y_, 10.0, 1e-12);
    EXPECT_NEAR(profile.shear_area_z_, 10.0, 1e-12);
}
