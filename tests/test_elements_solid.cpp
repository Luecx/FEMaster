/**
 * @file test_elements_solid.cpp
 * @brief Tests solid-element interpolation and compiled model result recovery.
 */

#include "../src/material/isotropic_elasticity.h"
#include "../src/model/model.h"
#include "../src/model/solid/c3d8.h"
#include "../src/model/solid/c3d8i.h"
#include "../src/section/section_solid.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace fem;

TEST(Elements_C3D8, ShapeFunctionsBasic) {
    model::C3D8 el(0, {0,1,2,3,4,5,6,7});
    const Precision pts[][3] = {{0,0,0},{-0.5,0.1,0.3},{0.7,-0.3,0.2}};
    for (auto &p : pts) {
        auto N = el.shape_function(p[0], p[1], p[2]);
        Precision s = 0;
        for (int i=0;i<8;++i) s += N(i,0);
        EXPECT_NEAR(s, 1.0, 1e-12);
    }

    const Precision corners[][3] = {
        {-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},
        {-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}
    };
    for (int n=0;n<8;++n){
        auto N = el.shape_function(corners[n][0], corners[n][1], corners[n][2]);
        for (int i=0;i<8;++i) EXPECT_NEAR(N(i,0), i==n?1.0:0.0, 1e-12);
    }
}

TEST(Elements_C3D8, ShapeDerivativeFD) {
    model::C3D8 el(0, {0,1,2,3,4,5,6,7});
    Precision r=0.2,s=-0.1,t=0.4, h=1e-6;
    auto dN = el.shape_derivative(r,s,t);

    auto N0 = el.shape_function(r,s,t);
    auto Nr = el.shape_function(r+h,s,t);
    auto Ns = el.shape_function(r,s+h,t);
    auto Nt = el.shape_function(r,s,t+h);
    for (int i=0;i<8;++i){
        EXPECT_NEAR((Nr(i,0)-N0(i,0))/h, dN(i,0), 1e-5);
        EXPECT_NEAR((Ns(i,0)-N0(i,0))/h, dN(i,1), 1e-5);
        EXPECT_NEAR((Nt(i,0)-N0(i,0))/h, dN(i,2), 1e-5);
    }
}

TEST(Elements_C3D8, TopAndBottomStressAreNotFlipped) {
    model::Model model;

    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_node(2, 1.0, 1.0, 0.0);
    model.set_node(3, 0.0, 1.0, 0.0);
    model.set_node(4, 0.0, 0.0, 1.0);
    model.set_node(5, 1.0, 0.0, 1.0);
    model.set_node(6, 1.0, 1.0, 1.0);
    model.set_node(7, 0.0, 1.0, 1.0);
    model.set_element<model::C3D8>(0, 0, 1, 2, 3, 4, 5, 6, 7);

    auto material = std::make_shared<material::Material>("MAT");
    material->set_elasticity<material::IsotropicElasticity>(1000.0, 0.3);
    model.add_material(material);

    const auto part = model._data->parts.get();
    ASSERT_NE(part, nullptr);

    auto section = std::make_shared<SolidSection>();
    section->material_ = material;
    section->region_   = part->elem_sets.get(SET_ELEM_ALL);
    model.add_section(section);
    model.compile();

    fem::model::Field displacement("U", fem::model::FieldDomain::NODE, 8, 6);
    displacement.set_zero();
    for (int node = 0; node < 8; ++node) {
        displacement(node, 0) = (node >= 4) ? 1.0 : 0.0;
    }

    auto [stress_top, stress_bot] = model.compute_stress_top_bot(displacement, false);

    ASSERT_EQ(stress_top.rows, stress_bot.rows);
    ASSERT_EQ(stress_top.components, stress_bot.components);
    for (Index i = 0; i < stress_top.rows; ++i) {
        for (Index j = 0; j < stress_top.components; ++j) {
            EXPECT_NEAR(stress_top(i, j), stress_bot(i, j), 1e-12);
        }
    }
}


TEST(Elements_C3D8I, AffinePatchMatchesC3D8OnDistortedHex) {
    model::Model model;

    const std::array<Vec3, 8> coords{
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.2, 0.0, 0.1),
        Vec3(1.1, 1.0, 0.0),
        Vec3(-0.1, 0.9, -0.1),
        Vec3(0.1, -0.1, 1.0),
        Vec3(1.1, 0.1, 1.2),
        Vec3(1.0, 1.1, 1.1),
        Vec3(0.0, 1.0, 0.9)
    };

    for (Index node = 0; node < 8; ++node) {
        model.set_node(node, coords[node].x(), coords[node].y(), coords[node].z());
    }

    model.set_element<model::C3D8>(0, 0, 1, 2, 3, 4, 5, 6, 7);
    model.set_element<model::C3D8I>(1, 0, 1, 2, 3, 4, 5, 6, 7);

    auto material = std::make_shared<material::Material>("MAT");
    material->set_elasticity<material::IsotropicElasticity>(210000.0, 0.3);
    model.add_material(material);

    auto section = std::make_shared<SolidSection>();
    section->material_ = material;
    section->region_ = model._data->parts.get()->elem_sets.get(SET_ELEM_ALL);
    model.add_section(section);

    model.compile();
    model.step_begin();

    auto* c3d8 = model._data->elements[0]->as<model::C3D8>();
    auto* c3d8i = model._data->elements[1]->as<model::C3D8I>();
    ASSERT_NE(c3d8, nullptr);
    ASSERT_NE(c3d8i, nullptr);
    EXPECT_EQ(c3d8i->type_name(), "C3D8I");
    EXPECT_EQ(c3d8i->num_ip(), 8);

    Precision c3d8_storage[24 * 24] {};
    Precision c3d8i_storage[24 * 24] {};
    const DynamicMatrix K = c3d8->stiffness(c3d8_storage);
    const DynamicMatrix KI = c3d8i->stiffness(c3d8i_storage);

    const Mat3 A = (Mat3() <<
        0.03,  0.01, -0.02,
       -0.01,  0.02,  0.015,
        0.005, -0.02, 0.01).finished();
    const Vec3 c(0.2, -0.1, 0.05);

    StaticVector<24> u = StaticVector<24>::Zero();
    for (Index node = 0; node < 8; ++node) {
        const Vec3 value = A * coords[node] + c;
        u.template segment<3>(3 * node) = value;
    }

    const DynamicVector force = K * u;
    const DynamicVector force_i = KI * u;
    EXPECT_LT((force_i - force).norm(), 1e-8 * (Precision(1) + force.norm()));

    // The condensed formulation must differ from the compatible C3D8 stiffness
    // away from the affine patch-test subspace.
    EXPECT_GT((KI - K).norm(), Precision(1e-6));

    model::Field displacement{
        "U", model::FieldDomain::NODE, 8, 6
    };
    displacement.set_zero();
    for (Index node = 0; node < 8; ++node) {
        displacement(node, 0) = u(3 * node + 0);
        displacement(node, 1) = u(3 * node + 1);
        displacement(node, 2) = u(3 * node + 2);
    }

    Precision kg_storage[24 * 24] {};
    const DynamicMatrix Kg = c3d8i->stiffness_geom(kg_storage, displacement);
    EXPECT_TRUE(Kg.allFinite());
    EXPECT_LT((Kg - Kg.transpose()).norm(), Precision(1e-10) * (Precision(1) + Kg.norm()));

    model.step_end();
}


TEST(Elements_C3D8I, NonlinearCondensedTangentMatchesFiniteDifference) {
    model::Model model;

    const std::array<Vec3, 8> coords{
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.2, 0.0, 0.1),
        Vec3(1.1, 1.0, 0.0),
        Vec3(-0.1, 0.9, -0.1),
        Vec3(0.1, -0.1, 1.0),
        Vec3(1.1, 0.1, 1.2),
        Vec3(1.0, 1.1, 1.1),
        Vec3(0.0, 1.0, 0.9)
    };

    for (Index node = 0; node < 8; ++node) {
        model.set_node(node, coords[node].x(), coords[node].y(), coords[node].z());
    }
    model.set_element<model::C3D8I>(0, 0, 1, 2, 3, 4, 5, 6, 7);

    auto material = std::make_shared<material::Material>("MAT");
    material->set_elasticity<material::IsotropicElasticity>(210000.0, 0.3);
    model.add_material(material);

    auto section = std::make_shared<SolidSection>();
    section->material_ = material;
    section->region_ = model._data->parts.get()->elem_sets.get(SET_ELEM_ALL);
    model.add_section(section);

    model.compile();
    model.step_begin();

    auto* element = model._data->elements[0]->as<model::C3D8I>();
    ASSERT_NE(element, nullptr);

    model::Field displacement{
        "U", model::FieldDomain::NODE, 8, 6
    };
    displacement.set_zero();

    for (Index node = 0; node < 8; ++node) {
        const Vec3& x = coords[node];
        displacement(node, 0) =  0.020 * x.x() * x.y();
        displacement(node, 1) = -0.015 * x.y() * x.z();
        displacement(node, 2) =  0.025 * x.x() * x.z();
    }

    model::NodeData internal{
        "INTERNAL_FORCES", model::FieldDomain::NODE, 8, 6
    };
    internal.set_zero();

    Precision storage[24 * 24] {};
    const DynamicMatrix tangent =
        element->stiffness_tangent(storage, internal, displacement);

    auto internal_force = [&](const model::Field& u) {
        model::NodeData force{
            "INTERNAL_FORCES", model::FieldDomain::NODE, 8, 6
        };
        force.set_zero();
        element->stiffness_tangent(nullptr, force, u);

        StaticVector<24> result = StaticVector<24>::Zero();
        for (Index node = 0; node < 8; ++node) {
            for (Dim dof = 0; dof < 3; ++dof) {
                result(3 * node + dof) = force(node, dof);
            }
        }
        return result;
    };

    constexpr Precision h = Precision(1e-7);
    const std::array<Index, 6> columns {0, 4, 8, 13, 17, 23};

    for (const Index column : columns) {
        model::Field plus = displacement;
        model::Field minus = displacement;

        const Index node = column / 3;
        const Index dof = column % 3;
        plus(node, dof) += h;
        minus(node, dof) -= h;

        const StaticVector<24> finite_difference =
            (internal_force(plus) - internal_force(minus)) / (Precision(2) * h);
        const DynamicVector tangent_column = tangent.col(column);

        EXPECT_LT(
            (finite_difference - tangent_column).norm(),
            Precision(2e-5) * (Precision(1) + finite_difference.norm())
        ) << "column=" << column;
    }

    // Verify objectivity under a superposed finite rigid rotation. The current
    // coordinates are rotated while the reference configuration remains fixed.
    const Precision angle = Precision(0.73);
    const Precision c     = std::cos(angle);
    const Precision s     = std::sin(angle);

    Mat3 Q;
    Q << c, -s, Precision(0),
         s,  c, Precision(0),
         Precision(0), Precision(0), Precision(1);

    model::Field rotated_displacement = displacement;
    for (Index node = 0; node < 8; ++node) {
        const Vec3 x_current = coords[node]
            + Vec3(
                displacement(node, 0),
                displacement(node, 1),
                displacement(node, 2)
            );
        const Vec3 u_rotated = Q * x_current - coords[node];

        rotated_displacement(node, 0) = u_rotated(0);
        rotated_displacement(node, 1) = u_rotated(1);
        rotated_displacement(node, 2) = u_rotated(2);
    }

    model::NodeData rotated_internal{
        "INTERNAL_FORCES", model::FieldDomain::NODE, 8, 6
    };
    rotated_internal.set_zero();

    Precision rotated_storage[24 * 24] {};
    const DynamicMatrix rotated_tangent =
        element->stiffness_tangent(
            rotated_storage,
            rotated_internal,
            rotated_displacement
        );

    StaticVector<24> base_force    = StaticVector<24>::Zero();
    StaticVector<24> rotated_force = StaticVector<24>::Zero();
    model::C3D8I::Matrix24 rotation_operator =
        model::C3D8I::Matrix24::Zero();

    for (Index node = 0; node < 8; ++node) {
        for (Dim dof = 0; dof < 3; ++dof) {
            base_force(3 * node + dof) = internal(node, dof);
            rotated_force(3 * node + dof) = rotated_internal(node, dof);
        }

        rotation_operator.block<3, 3>(3 * node, 3 * node) = Q;
    }

    const StaticVector<24> expected_rotated_force =
        rotation_operator * base_force;
    const model::C3D8I::Matrix24 expected_rotated_tangent =
        rotation_operator * tangent * rotation_operator.transpose();

    EXPECT_LT(
        (rotated_force - expected_rotated_force).norm(),
        Precision(1e-9) * (Precision(1) + base_force.norm())
    );
    EXPECT_LT(
        (rotated_tangent - expected_rotated_tangent).norm(),
        Precision(1e-8) * (Precision(1) + tangent.norm())
    );

    // Model-level output uses the thermal-aware virtual overload even when no
    // thermal field is supplied; this call therefore also checks C3D8I recovery
    // dispatch through StructuralElement.
    EXPECT_NO_THROW(model.compute_stress_nodal(displacement, true));

    model.step_end();
}
