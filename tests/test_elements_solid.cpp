/**
 * @file test_elements_solid.cpp
 * @brief Tests compatible and enhanced solid-element mechanics.
 *
 * The checks cover C3D8 interpolation and C3D8I affine consistency, condensed
 * finite-strain tangents, rigid-motion objectivity and the linear rigid-body
 * nullspace. Plastic loading exercises committed/trial history through the
 * common nonlinear state manager and verifies that auxiliary element paths
 * remain state-neutral.
 */

#include "../src/material/isotropic_elasticity.h"
#include "../src/material/isotropic_j2_elasticity.h"
#include "../src/material/neo_hooke_elasticity.h"
#include "../src/loadcase/tools/nonlinear_state_manager.h"
#include "../src/model/model.h"
#include "../src/model/solid/c3d8.h"
#include "../src/model/solid/c3d8i.h"
#include "../src/section/section_solid.h"

#include <gtest/gtest.h>

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

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
    const DynamicMatrix K = c3d8->evaluate(c3d8_storage, nullptr, nullptr, nullptr, nullptr, nullptr, false);
    const DynamicMatrix KI = c3d8i->evaluate(c3d8i_storage, nullptr, nullptr, nullptr, nullptr, nullptr, false);

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
    const DynamicMatrix Kg = c3d8i->evaluate(nullptr, kg_storage, nullptr, &displacement, nullptr, nullptr, false);
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
        element->evaluate(storage, nullptr, &internal, &displacement, &displacement, nullptr, true);

    auto internal_force = [&](const model::Field& u) {
        model::NodeData force{
            "INTERNAL_FORCES", model::FieldDomain::NODE, 8, 6
        };
        force.set_zero();
        element->evaluate(nullptr, nullptr, &force, &u, &u, nullptr, true);

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
        element->evaluate(
            rotated_storage,
            nullptr,
            &rotated_internal,
            &rotated_displacement,
            &rotated_displacement,
            nullptr,
            true
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

/**
 * Checks every condensed tangent column against centered differences of the
 * stationary internal force for two hyperelastic laws and distorted geometry.
 *
 * The perturbed evaluations solve the local enhanced equations independently.
 * Thus the comparison includes the implicit enhanced-state derivative rather
 * than holding the local parameters fixed. A finite spatial rotation and
 * translation additionally exercise force and tangent objectivity.
 */
TEST(Elements_C3D8I, CompleteHyperelasticTangentsAndRigidMotions) {
    Precision maximum_tangent_error         = Precision(0);
    Precision maximum_rotation_force_error  = Precision(0);
    Precision maximum_rotation_tangent_error = Precision(0);

    const std::array<Vec3, 8> regular {
        Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(1, 1, 0), Vec3(0, 1, 0),
        Vec3(0, 0, 1), Vec3(1, 0, 1), Vec3(1, 1, 1), Vec3(0, 1, 1)
    };
    const std::array<Vec3, 8> distorted {
        Vec3(0, 0, 0), Vec3(1.2, 0, 0.1), Vec3(1.1, 1, 0), Vec3(-0.1, 0.9, -0.1),
        Vec3(0.1, -0.1, 1), Vec3(1.1, 0.1, 1.2), Vec3(1, 1.1, 1.1), Vec3(0, 1, 0.9)
    };

    // Exercise compressible and nearly incompressible tangents for both laws
    for (const auto& coords : {regular, distorted}) {
        for (const bool neo_hooke : {false, true}) {
            for (const Precision poisson : {Precision(0.3), Precision(0.4999)}) {
                model::Model model;
                for (Index node = 0; node < 8; ++node) {
                    model.set_node(node, coords[node].x(), coords[node].y(), coords[node].z());
                }
                model.set_element<model::C3D8I>(0, 0, 1, 2, 3, 4, 5, 6, 7);

                auto material = std::make_shared<material::Material>("MAT");
                if (neo_hooke) {
                    const Precision shear = Precision(210000) / (Precision(2) * (Precision(1) + poisson));
                    const Precision bulk  = Precision(210000) / (Precision(3) * (Precision(1) - Precision(2) * poisson));
                    material->set_elasticity<material::NeoHookeElasticity>(shear / Precision(2), Precision(2) / bulk);
                } else {
                    material->set_elasticity<material::IsotropicElasticity>(210000, poisson);
                }
                model.add_material(material);

                auto section = std::make_shared<SolidSection>();
                section->material_ = material;
                section->region_   = model._data->parts.get()->elem_sets.get(SET_ELEM_ALL);
                model.add_section(section);
                model.compile();
                model.step_begin();

                auto* element = model._data->elements[0]->as<model::C3D8I>();
                ASSERT_NE(element, nullptr);

                // Reconstruct the force through the same public residual-only path
                auto internal_force = [&](const model::Field& displacement) {
                    model::NodeData forces("F", model::FieldDomain::NODE, 8, 6);
                    forces.set_zero();
                    element->evaluate(nullptr, nullptr, &forces, &displacement, &displacement, nullptr, true);

                    StaticVector<24> result;
                    for (Index node = 0; node < 8; ++node) {
                        for (Dim component = 0; component < 3; ++component) {
                            result(3 * node + component) = forces(node, component);
                        }
                    }
                    return result;
                };

                for (const Precision amplitude : {Precision(-0.10), Precision(0.02), Precision(0.15)}) {
                    SCOPED_TRACE(::testing::Message() << "neo_hooke=" << neo_hooke
                        << ", nu=" << poisson << ", amplitude=" << amplitude
                        << ", distorted=" << (coords[1].x() != Precision(1)));

                    // Combine finite affine deformation with non-affine bending
                    model::Field displacement("U", model::FieldDomain::NODE, 8, 6);
                    displacement.set_zero();
                    for (Index node = 0; node < 8; ++node) {
                        const Vec3& x = coords[node];
                        displacement(node, 0) = amplitude * (Precision(0.3) * x.x() + x.x() * x.y());
                        displacement(node, 1) = amplitude * (Precision(-0.2) * x.y() - x.y() * x.z());
                        displacement(node, 2) = amplitude * (Precision(0.1) * x.z() + x.x() * x.z());
                    }

                    model::NodeData forces("F", model::FieldDomain::NODE, 8, 6);
                    forces.set_zero();
                    Precision storage[24 * 24] {};
                    const DynamicMatrix tangent = element->evaluate(storage, nullptr, &forces, &displacement, &displacement, nullptr, true);
                    EXPECT_LT((tangent - tangent.transpose()).norm(), Precision(1e-10) * tangent.norm());

                    // Difference every independent translational degree of freedom
                    DynamicMatrix difference(24, 24);
                    constexpr Precision h = Precision(2e-7);
                    for (Index column = 0; column < 24; ++column) {
                        model::Field plus  = displacement;
                        model::Field minus = displacement;
                        plus (column / 3, column % 3) += h;
                        minus(column / 3, column % 3) -= h;
                        difference.col(column) = (internal_force(plus) - internal_force(minus)) / (Precision(2) * h);
                    }
                    const Precision derivative_error = (difference - tangent).norm() / tangent.norm();
                    EXPECT_LT(derivative_error, Precision(2e-6));
                    maximum_tangent_error = std::max(maximum_tangent_error, derivative_error);

                    // Superpose an arbitrary finite rotation and translation
                    const Mat3 Q = Eigen::AngleAxis<Precision>(Precision(1.2), Vec3(1, -2, 3).normalized()).toRotationMatrix();
                    const Vec3 translation(0.7, -0.4, 0.3);
                    model::Field rotated = displacement;
                    StaticMatrix<24, 24> rotation = StaticMatrix<24, 24>::Zero();
                    for (Index node = 0; node < 8; ++node) {
                        const Vec3 current = coords[node] + Vec3(displacement(node, 0), displacement(node, 1), displacement(node, 2));
                        const Vec3 value   = Q * current + translation - coords[node];
                        for (Dim component = 0; component < 3; ++component) {
                            rotated(node, component) = value(component);
                        }
                        rotation.block<3, 3>(3 * node, 3 * node) = Q;
                    }
                    forces.set_zero();
                    const DynamicMatrix rotated_tangent = element->evaluate(storage, nullptr, &forces, &rotated, &rotated, nullptr, true);
                    const auto base_force    = internal_force(displacement);
                    const auto rotated_force = internal_force(rotated);
                    const Precision force_error = (rotated_force - rotation * base_force).norm()
                        / (Precision(1) + base_force.norm());
                    const Precision rotation_error = (rotated_tangent - rotation * tangent * rotation.transpose()).norm()
                        / tangent.norm();
                    EXPECT_LT(force_error, Precision(1e-8));
                    EXPECT_LT(rotation_error, Precision(1e-9));
                    maximum_rotation_force_error   = std::max(maximum_rotation_force_error, force_error);
                    maximum_rotation_tangent_error = std::max(maximum_rotation_tangent_error, rotation_error);
                }
                model.step_end();
            }
        }
    }

    // Retain the worst errors from every material, geometry and deformation case
    RecordProperty("maximum_tangent_relative_error", (::testing::Message() << maximum_tangent_error).GetString());
    RecordProperty("maximum_rotation_force_relative_error", (::testing::Message() << maximum_rotation_force_error).GetString());
    RecordProperty("maximum_rotation_tangent_relative_error", (::testing::Message() << maximum_rotation_tangent_error).GetString());
}

/**
 * Checks the undeformed tangent, rigid-body nullspace and linear stiffness for
 * an element spanning a range of volumetric-to-shear stiffness ratios.
 */
TEST(Elements_C3D8I, LinearLimitAndSixRigidBodyModes) {
    for (const Precision poisson : {Precision(0), Precision(0.3), Precision(0.4999)}) {
        SCOPED_TRACE(poisson);
        model::Model model;
        const std::array<Vec3, 8> coords {
            Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(1, 1, 0), Vec3(0, 1, 0),
            Vec3(0, 0, 1), Vec3(1, 0, 1), Vec3(1, 1, 1), Vec3(0, 1, 1)
        };
        for (Index node = 0; node < 8; ++node) {
            model.set_node(node, coords[node].x(), coords[node].y(), coords[node].z());
        }
        model.set_element<model::C3D8I>(0, 0, 1, 2, 3, 4, 5, 6, 7);
        auto material = std::make_shared<material::Material>("MAT");
        material->set_elasticity<material::IsotropicElasticity>(210000, poisson);
        model.add_material(material);
        auto section = std::make_shared<SolidSection>();
        section->material_ = material;
        section->region_   = model._data->parts.get()->elem_sets.get(SET_ELEM_ALL);
        model.add_section(section);
        model.compile();
        model.step_begin();

        // Compare linear condensation with the undeformed finite-strain tangent
        auto* element = model._data->elements[0]->as<model::C3D8I>();
        model::Field displacement("U", model::FieldDomain::NODE, 8, 6);
        model::NodeData forces("F", model::FieldDomain::NODE, 8, 6);
        displacement.set_zero();
        forces.set_zero();
        Precision storage[24 * 24] {};
        const DynamicMatrix linear    = element->evaluate(storage, nullptr, nullptr, nullptr, nullptr, nullptr, false);
        const DynamicMatrix nonlinear = element->evaluate(storage, nullptr, &forces, &displacement, &displacement, nullptr, true);
        EXPECT_LT((linear - nonlinear).norm(), Precision(1e-12) * linear.norm());
        EXPECT_LT((linear - linear.transpose()).norm(), Precision(1e-14) * linear.norm());

        // A stable eight-node continuum has six rigid modes and eighteen positive modes
        Eigen::SelfAdjointEigenSolver<DynamicMatrix> eigenvalues(linear);
        ASSERT_EQ(eigenvalues.info(), Eigen::Success);
        for (Index mode = 0; mode < 6; ++mode) {
            EXPECT_LT(std::abs(eigenvalues.eigenvalues()(mode)), Precision(1e-10) * linear.norm());
        }
        EXPECT_GT(eigenvalues.eigenvalues()(6), Precision(1));
        model.step_end();
    }
}

/**
 * Exercises plastic loading, commitment, unloading and state-neutral auxiliary
 * element evaluations. Only a physical trial force/tangent evaluation may
 * overwrite material_state_new; material_state_old remains immutable.
 */
TEST(Elements_C3D8I, PlasticHistoryAndAuxiliaryStateNeutrality) {
    model::Model model;
    const std::array<Vec3, 8> coords {
        Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(1, 1, 0), Vec3(0, 1, 0),
        Vec3(0, 0, 1), Vec3(1, 0, 1), Vec3(1, 1, 1), Vec3(0, 1, 1)
    };
    for (Index node = 0; node < 8; ++node) {
        model.set_node(node, coords[node].x(), coords[node].y(), coords[node].z());
    }
    model.set_element<model::C3D8I>(0, 0, 1, 2, 3, 4, 5, 6, 7);
    auto material = std::make_shared<material::Material>("MAT");
    material->set_elasticity<material::IsotropicJ2Elasticity>(210000, 0.3);
    auto* j2 = dynamic_cast<material::IsotropicJ2Elasticity*>(material->elasticity().get());
    j2->add_yield_point(280, 0);
    j2->add_yield_point(350, 0.1);
    model.add_material(material);
    auto section = std::make_shared<SolidSection>();
    section->material_ = material;
    section->region_   = model._data->parts.get()->elem_sets.get(SET_ELEM_ALL);
    model.add_section(section);
    model.compile();
    model.step_begin();
    loadcase::tools::NonlinearStateManager state(model);

    auto* element = model._data->elements[0]->as<model::C3D8I>();

    // Commit two tensile increments and then unload from the second committed state
    for (const Precision stretch : {Precision(0.004), Precision(0.008), Precision(0.007)}) {
        SCOPED_TRACE(stretch);
        model::Field displacement("U", model::FieldDomain::NODE, 8, 6);
        displacement.set_zero();
        for (Index node = 0; node < 8; ++node) {
            displacement(node, 0) = stretch * coords[node].x();
            displacement(node, 1) = -Precision(0.3) * stretch * coords[node].y();
            displacement(node, 2) = -Precision(0.3) * stretch * coords[node].z();
            displacement(node, 0) += Precision(0.15) * stretch * coords[node].x() * coords[node].y();
            displacement(node, 2) += Precision(0.10) * stretch * coords[node].x() * coords[node].z();
        }

        // Physical trial evaluations must be deterministic from the committed history
        const model::Field committed = *model._data->material_state_old;
        model::NodeData forces("F", model::FieldDomain::NODE, 8, 6);
        forces.set_zero();
        Precision storage[24 * 24] {};
        const DynamicMatrix tangent = element->evaluate(storage, nullptr, &forces, &displacement, &displacement, nullptr, true);
        const model::Field trial = *model._data->material_state_new;
        const model::Field original_force = forces;
        forces.set_zero();
        element->evaluate(nullptr, nullptr, &forces, &displacement, &displacement, nullptr, true);
        EXPECT_TRUE(tangent.allFinite());

        // Auxiliary paths may reconstruct local parameters but must not commit history
        element->evaluate(storage, nullptr, nullptr, nullptr, nullptr, nullptr, false);
        element->evaluate(nullptr, storage, nullptr, &displacement, nullptr, nullptr, false);
        EXPECT_NO_THROW(model.compute_stress_nodal(displacement, true));
        EXPECT_NO_THROW(model.compute_stress_nodal(displacement, false));
        for (Index row = 0; row < trial.rows; ++row) {
            for (Dim component = 0; component < trial.components; ++component) {
                EXPECT_DOUBLE_EQ((*model._data->material_state_old)(row, component), committed(row, component));
                EXPECT_NEAR((*model._data->material_state_new)(row, component), trial(row, component), Precision(1e-12));
            }
            EXPECT_GE(trial(row, 6), committed(row, 6) - Precision(1e-12));
        }
        for (Index node = 0; node < 8; ++node) {
            for (Dim component = 0; component < 3; ++component) {
                EXPECT_NEAR(forces(node, component), original_force(node, component), Precision(1e-9));
            }
        }
        EXPECT_GT(trial(0, 6), Precision(0));
        state.commit_material_state();
    }
    model.step_end();
}
