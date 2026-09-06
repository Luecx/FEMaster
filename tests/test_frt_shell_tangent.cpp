/**
 * @file test_frt_shell_tangent.cpp
 * @brief Smoke-tests finite-rotation shell tangents for all supported topologies.
 */

#include "../src/material/isotropic_elasticity.h"
#include "../src/model/element/element_structural.h"
#include "../src/model/model.h"
#include "../src/model/shell/frt_shell_s3.h"
#include "../src/model/shell/frt_shell_s4.h"
#include "../src/model/shell/frt_shell_s6.h"
#include "../src/model/shell/frt_shell_s8.h"
#include "../src/section/section_shell_integrated.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

namespace {

using namespace fem;

template<typename Element, std::size_t... I>
void add_element(model::Model& model, std::index_sequence<I...>) {
    model.set_element<Element>(0, static_cast<ID>(I)...);
}

template<typename Element, std::size_t N>
void expect_frt_tangent_smoke(const std::array<Vec3, N>& nodes) {
    model::Model model;

    for (std::size_t i = 0; i < N; ++i)
        model.set_node(
            static_cast<ID>(i),
            nodes[i](0), nodes[i](1), nodes[i](2));

    add_element<Element>(model, std::make_index_sequence<N>{});

    auto material = std::make_shared<material::Material>("MAT");
    material->set_elasticity<material::IsotropicElasticity>(1000.0, 0.3);
    model.add_material(material);

    model.add_section(std::make_shared<IntegratedShellSection>(
        material,
        model._data->parts.get()->elem_sets.get(SET_ELEM_ALL),
        Precision(0.1),
        nullptr
    ));

    model.compile();
    model.step_begin();

    auto* structural = model._data->elements[0]->as<model::StructuralElement>();
    ASSERT_NE(structural, nullptr);

    model::Field displacement{
        "U", model::FieldDomain::NODE, static_cast<Index>(N), 6};
    displacement.set_zero();

    model::NodeData nodal_forces{
        "INTERNAL_FORCE", model::FieldDomain::NODE, static_cast<Index>(N), 6};
    nodal_forces.set_zero();

    constexpr Index ndofs = static_cast<Index>(6 * N);
    std::vector<Precision> storage(static_cast<std::size_t>(ndofs * ndofs));
    const DynamicMatrix tangent = structural->stiffness_tangent(
        storage.data(), nodal_forces, displacement);

    ASSERT_EQ(tangent.rows(), ndofs);
    ASSERT_EQ(tangent.cols(), ndofs);
    EXPECT_TRUE(tangent.allFinite());

    const Precision scale = std::max(Precision(1), tangent.norm());
    EXPECT_LT((tangent - tangent.transpose()).norm() / scale, Precision(1e-10));

    for (Index node = 0; node < static_cast<Index>(N); ++node)
        for (Index component = 0; component < 6; ++component)
            EXPECT_TRUE(std::isfinite(nodal_forces(node, component)));

    model.step_end();
}

} // namespace

TEST(FRTShellTangentSmoke, S3) {
    const std::array<Vec3, 3> nodes{{
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.05),
        Vec3(0.0, 1.0, -0.02)
    }};
    expect_frt_tangent_smoke<model::FRTShellS3>(nodes);
}

TEST(FRTShellTangentSmoke, S4) {
    const std::array<Vec3, 4> nodes{{
        Vec3(0.0, 0.0, 0.00),
        Vec3(1.0, 0.0, 0.05),
        Vec3(1.0, 1.0, 0.12),
        Vec3(0.0, 1.0, -0.03)
    }};
    expect_frt_tangent_smoke<model::FRTShellS4>(nodes);
}

TEST(FRTShellTangentSmoke, S6) {
    const std::array<Vec3, 6> nodes{{
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(0.5, 0.0, 0.0),
        Vec3(0.5, 0.5, 0.0),
        Vec3(0.0, 0.5, 0.0)
    }};
    expect_frt_tangent_smoke<model::FRTShellS6>(nodes);
}

TEST(FRTShellTangentSmoke, S8) {
    const std::array<Vec3, 8> nodes{{
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0),
        Vec3(1.0, 1.0, 0.0),
        Vec3(0.0, 1.0, 0.0),
        Vec3(0.5, 0.0, 0.0),
        Vec3(1.0, 0.5, 0.0),
        Vec3(0.5, 1.0, 0.0),
        Vec3(0.0, 0.5, 0.0)
    }};
    expect_frt_tangent_smoke<model::FRTShellS8>(nodes);
}
