/**
 * @file test_shell_thermal.cpp
 * @brief Uniform-through-thickness thermal equivalent loads for MITC/FRT shells.
 */

#include "../src/bc/neumann/load_t.h"
#include "../src/material/isotropic_elasticity.h"
#include "../src/material/isotropic_j2_elasticity.h"
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
#include <tuple>
#include <utility>

namespace {

template<typename Shell, std::size_t... I>
void add_shell(fem::model::Model& model, std::index_sequence<I...>) {
    model.set_element<Shell>(0, static_cast<fem::ID>(I)...);
}

/**
 * Checks uniform shell thermal loading against free isotropic dilation.
 *
 * Loads must equal K*u for the manufactured expansion, and integrated-section
 * recovery must produce zero stress, zero resultants and zero prestress
 * stiffness. Plane patches also check the analytical restrained thermal stress.
 */
template<typename Shell, std::size_t N>
void check_uniform_thermal_load(
    const std::array<fem::Vec3, N>& coords, bool use_abd, fem::Precision free_strain = 0.2
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
        (*temperature)(i, 0) = 20.0 + free_strain / 0.01;
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
    const DynamicMatrix K = element->evaluate(matrix_storage, nullptr, nullptr, nullptr, nullptr, nullptr, false);
    StaticVector<6 * N> free_expansion = StaticVector<6 * N>::Zero();
    for (Index node = 0; node < static_cast<Index>(N); ++node) {
        free_expansion(6 * node + 0) = free_strain * coords[node].x();
        free_expansion(6 * node + 1) = free_strain * coords[node].y();
        free_expansion(6 * node + 2) = free_strain * coords[node].z();
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

    if (!use_abd) {
        model::Field thermal_free_strain{
            "THERMAL_FREE_STRAIN",
            model::FieldDomain::ELEMENT_NODAL,
            model._data->field_rows(model::FieldDomain::ELEMENT_NODAL),
            1
        };
        thermal_free_strain.set_zero();
        load.apply_thermal_free_strain(*model._data, thermal_free_strain);

        model::Field displacement{
            "DISPLACEMENT", model::FieldDomain::NODE, static_cast<Index>(N), 6
        };
        displacement.set_zero();
        for (Index node = 0; node < static_cast<Index>(N); ++node) {
            displacement(node, 0) = free_strain * coords[static_cast<std::size_t>(node)].x();
            displacement(node, 1) = free_strain * coords[static_cast<std::size_t>(node)].y();
            displacement(node, 2) = free_strain * coords[static_cast<std::size_t>(node)].z();
        }

        auto stress_strain =
            model.compute_stress_nodal(displacement, nullptr, &thermal_free_strain);
        const auto& stress = std::get<0>(stress_strain);
        const auto& strain = std::get<1>(stress_strain);

        bool planar_xy = true;
        for (Index node = 1; node < static_cast<Index>(N); ++node) {
            planar_xy = planar_xy
                && std::abs(coords[static_cast<std::size_t>(node)].z() - coords[0].z()) < 1e-12;
        }

        for (Index node = 0; node < static_cast<Index>(N); ++node) {
            if (planar_xy) {
                EXPECT_NEAR(strain(node, 0), free_strain, 1e-9);
                EXPECT_NEAR(strain(node, 1), free_strain, 1e-9);
            }
            for (Index component = 0; component < 6; ++component) {
                EXPECT_NEAR(stress(node, component), 0.0, 1e-8)
                    << "node=" << node << ", component=" << component;
            }
        }

        const auto resultants =
            model.compute_shell_resultants(displacement, &thermal_free_strain);
        for (Index node = 0; node < static_cast<Index>(N); ++node) {
            for (Index component = 0; component < 8; ++component) {
                EXPECT_NEAR(resultants(node, component), 0.0, 1e-8)
                    << "node=" << node << ", resultant=" << component;
            }
        }

        Precision kg_storage[6 * N * 6 * N] {};
        const DynamicMatrix Kg_free =
            element->evaluate(nullptr, kg_storage, nullptr, &displacement, nullptr, &thermal_free_strain, false);
        EXPECT_LT(Kg_free.norm(), 1e-8);

        displacement.set_zero();

        // Fully restrained planar expansion has equal biaxial compression:
        // sigma_x = sigma_y = -E * alpha * DeltaT / (1 - nu).
        if (planar_xy) {
            const auto restrained =
                model.compute_stress_nodal(displacement, nullptr, &thermal_free_strain);
            const auto& restrained_stress = std::get<0>(restrained);
            for (Index node = 0; node < static_cast<Index>(N); ++node) {
                EXPECT_NEAR(restrained_stress(node, 0), -1000.0 * free_strain / 0.75, 1e-8);
                EXPECT_NEAR(restrained_stress(node, 1), -1000.0 * free_strain / 0.75, 1e-8);
                for (Index component = 2; component < 6; ++component) {
                    EXPECT_NEAR(restrained_stress(node, component), 0.0, 1e-8);
                }
            }
        }

        const DynamicMatrix Kg_restrained =
            element->evaluate(nullptr, kg_storage, nullptr, &displacement, nullptr, &thermal_free_strain, false);
        EXPECT_GT(Kg_restrained.norm(), 1e-6);
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

TEST(ShellPlasticity, PeeqUsesThicknessMaximumBeforeExtrapolation) {
    using namespace fem;

    model::Model model;
    model.set_node(0, 0.0, 0.0, 0.0);
    model.set_node(1, 1.0, 0.0, 0.0);
    model.set_node(2, 1.0, 1.0, 0.0);
    model.set_node(3, 0.0, 1.0, 0.0);
    model.set_element<model::FRTShellS4>(0, 0, 1, 2, 3);

    auto material = std::make_shared<material::Material>("MAT");
    material->set_elasticity<material::IsotropicJ2Elasticity>(210000.0, 0.3);
    model.add_material(material);

    const auto region = model._data->parts.get()->elem_sets.get(SET_ELEM_ALL);
    model.add_section(std::make_shared<IntegratedShellSection>(
        material, region, Precision(0.2), nullptr
    ));
    model.compile();

    auto* shell = model._data->elements[0]->as<model::FRTShellS4>();
    ASSERT_NE(shell, nullptr);

    // Expand committed state to the J2 layout and initialize every shell material point
    const Index material_points = model._data->field_rows(model::FieldDomain::ELEMENT_MP);
    model._data->material_state_old = std::make_shared<model::Field>(
        "MATERIAL_STATE_OLD", model::FieldDomain::ELEMENT_MP, material_points, 7
    );
    model.initialize_material_state(*model._data->material_state_old);

    const RowMatrix ip_rst = shell->stress_strain_ip_rst();

    // Prescribe a bilinear in-plane PEEQ field. At each in-plane IP exactly one
    // of the five thickness points carries the target maximum; all others are lower.
    for (Index ip = 0; ip < static_cast<Index>(ip_rst.rows()); ++ip) {
        const Precision r = ip_rst(ip, 0);
        const Precision s = ip_rst(ip, 1);
        const Precision expected = Precision(0.2) + Precision(0.03) * r
                                 + Precision(0.04) * s + Precision(0.02) * r * s;

        for (Index mp = 0; mp < shell->num_mp_per_ip(); ++mp) {
            (*model._data->material_state_old)(shell->mp_index(ip, mp), 6) =
                expected - Precision(0.01) * static_cast<Precision>(mp + 1);
        }

        const Index maximum_mp = (ip + 2) % shell->num_mp_per_ip();
        (*model._data->material_state_old)(shell->mp_index(ip, maximum_mp), 6) = expected;
    }

    // The shell first reduces through thickness and then extrapolates the IP
    // maxima to nodes. A bilinear source field must therefore be reproduced exactly.
    const model::Field peeq = model.compute_peeq_nodal();
    const RowMatrix nodal_rst = shell->stress_strain_nodal_rst();

    for (Index node = 0; node < static_cast<Index>(nodal_rst.rows()); ++node) {
        const Precision r = nodal_rst(node, 0);
        const Precision s = nodal_rst(node, 1);
        const Precision expected = Precision(0.2) + Precision(0.03) * r
                                 + Precision(0.04) * s + Precision(0.02) * r * s;
        EXPECT_NEAR(peeq(node, 0), expected, Precision(1e-12));
    }
}

TEST(ShellThermal, S3Integrated) {
    check_uniform_thermal_load<fem::model::FRTShellS3, 3>(
        {fem::Vec3(0, 0, 0), fem::Vec3(1, 0, 0), fem::Vec3(0, 1, 0)}, false
    );
}

TEST(ShellThermal, S4Integrated) {
    check_uniform_thermal_load<fem::model::FRTShellS4, 4>(
        {fem::Vec3(0, 0, 0), fem::Vec3(1, 0, 0),
         fem::Vec3(1, 1, 0), fem::Vec3(0, 1, 0)}, false
    );
}

TEST(ShellThermal, S6Integrated) {
    check_uniform_thermal_load<fem::model::FRTShellS6, 6>(
        {fem::Vec3(0, 0, 0), fem::Vec3(1, 0, 0), fem::Vec3(0, 1, 0),
         fem::Vec3(0.5, 0, 0), fem::Vec3(0.5, 0.5, 0),
         fem::Vec3(0, 0.5, 0)}, false
    );
}

TEST(ShellThermal, S8Integrated) {
    check_uniform_thermal_load<fem::model::FRTShellS8, 8>(
        {fem::Vec3(0, 0, 0), fem::Vec3(1, 0, 0),
         fem::Vec3(1, 1, 0), fem::Vec3(0, 1, 0),
         fem::Vec3(0.5, 0, 0), fem::Vec3(1, 0.5, 0),
         fem::Vec3(0.5, 1, 0), fem::Vec3(0, 0.5, 0)}, false
    );
}

TEST(ShellThermal, S4ABDWithMembraneBendingCoupling) {
    check_uniform_thermal_load<fem::model::FRTShellS4, 4>(
        {fem::Vec3(0, 0, 0), fem::Vec3(1, 0, 0),
         fem::Vec3(1, 1, 0), fem::Vec3(0, 1, 0)}, true
    );
}

TEST(ShellThermal, S8CurvedIntegrated) {
    constexpr fem::Precision radius = 2.0;
    constexpr fem::Precision angle = 0.2;

    const auto p = [=](fem::Precision theta, fem::Precision y) {
        return fem::Vec3(
            radius * std::sin(theta),
            y,
            radius * std::cos(theta)
        );
    };

    check_uniform_thermal_load<fem::model::FRTShellS8, 8>(
        {p(-angle, 0.0), p( angle, 0.0), p( angle, 1.0), p(-angle, 1.0),
         p( 0.0,   0.0), p( angle, 0.5), p( 0.0,   1.0), p(-angle, 0.5)},
        false
    );
}

/**
 * Checks uniform dilation across curvature, thermal sign and global orientation.
 *
 * Curved MITC8 strains and thermal initial strains must use the same discrete
 * mapping for both integrated and membrane-bending-coupled ABD sections. A
 * rigidly rotated and translated patch exercises the physical tangent bases.
 */
TEST(ShellThermal, S8CurvedExpansionSweep) {
    using namespace fem;

    // Sample gentle and strong curvature and both heating and cooling
    for (const Precision angle : {Precision(0.02), Precision(0.2), Precision(0.6)}) {
        for (const Precision radius : {Precision(0.7), Precision(2), Precision(10)}) {
            for (const Precision free_strain : {Precision(0.001), Precision(-0.002)}) {
                for (const bool use_abd : {false, true}) {
                    SCOPED_TRACE(::testing::Message()
                        << "angle=" << angle << ", radius=" << radius
                        << ", free_strain=" << free_strain << ", ABD=" << use_abd);

                    // Rotate the cylinder about two global axes and move its origin
                    const Precision c = std::cos(Precision(0.4));
                    const Precision s = std::sin(Precision(0.4));
                    Mat3 rotation;
                    rotation << c, -s*c, s*s,
                                s,  c*c, -c*s,
                                0,    s,    c;
                    const auto p = [&](Precision theta, Precision y) -> Vec3 {
                        const Vec3 local(radius * std::sin(theta), y, radius * std::cos(theta));
                        return rotation * local + Vec3(3, -2, 1);
                    };

                    check_uniform_thermal_load<model::FRTShellS8, 8>(
                        {p(-angle, 0), p(angle, 0), p(angle, 1), p(-angle, 1),
                         p(0, 0), p(angle, 0.5), p(0, 1), p(-angle, 0.5)},
                        use_abd, free_strain
                    );
                }
            }
        }
    }
}

/**
 * Checks a compatible nonuniform thermal expansion on a planar S8 patch.
 *
 * For epsilon_th(x) = a + b*x the quadratic displacement
 * u_x = a*x + b*(x*x-y*y)/2, u_y = a*y + b*x*y has exactly the prescribed
 * isotropic strain and zero engineering shear. Its infinitesimal drilling
 * rotation is b*y. This verifies temperature interpolation independently of
 * uniform dilation and of the thermal implementation's MITC mapping.
 */
TEST(ShellThermal, S8PlanarLinearTemperatureGradient) {
    using namespace fem;

    // Build a planar serendipity patch with a known quadratic free expansion
    const std::array<Vec3, 8> coords{
        Vec3(-1,-1,0), Vec3(1,-1,0), Vec3(1,1,0), Vec3(-1,1,0),
        Vec3(0,-1,0), Vec3(1,0,0), Vec3(0,1,0), Vec3(-1,0,0)
    };
    model::Model model;
    for (Index node = 0; node < 8; ++node) {
        model.set_node(node, coords[node].x(), coords[node].y(), coords[node].z());
    }
    model.set_element<model::FRTShellS8>(0, 0,1,2,3,4,5,6,7);
    auto material = std::make_shared<material::Material>("MAT");
    material->set_elasticity<material::IsotropicElasticity>(1000.0, 0.25);
    material->set_thermal_expansion(0.01);
    const auto region = model._data->parts.get()->elem_sets.get(SET_ELEM_ALL);
    model.add_material(material);
    model.add_section(std::make_shared<IntegratedShellSection>(material, region, 0.2, nullptr));
    model.compile();
    model.step_begin();

    // Match the linear temperature field with its analytical compatible motion
    constexpr Precision a = 0.001;
    constexpr Precision b = 0.0003;
    auto temperature = std::make_shared<model::Field>("TEMP", model::FieldDomain::NODE, 8, 1);
    model::Field displacement{"DISPLACEMENT", model::FieldDomain::NODE, 8, 6};
    displacement.set_zero();
    StaticVector<48> motion = StaticVector<48>::Zero();
    for (Index node = 0; node < 8; ++node) {
        const Precision x = coords[node].x();
        const Precision y = coords[node].y();
        (*temperature)(node, 0) = 20.0 + (a + b*x) / 0.01;
        displacement(node, 0) = a*x + b*(x*x-y*y)/2;
        displacement(node, 1) = a*y + b*x*y;
        displacement(node, 5) = b*y;
        for (Index dof = 0; dof < 6; ++dof) {
            motion(6*node+dof) = displacement(node, dof);
        }
    }

    bc::TLoad load;
    load.temp_field_ = temperature;
    load.ref_temp_   = 20.0;
    model::Field rhs{"RHS", model::FieldDomain::NODE, 8, 6};
    rhs.set_zero();
    load.apply(*model._data, rhs, 0.0);

    // Verify consistent nodal loads, including the independent drilling motion
    auto* element = model._data->elements[0]->as<model::FRTShellS8>();
    ASSERT_NE(element, nullptr);
    Precision storage[48*48]{};
    const DynamicMatrix K = element->evaluate(storage, nullptr, nullptr, nullptr, nullptr, nullptr, false);
    const StaticVector<48> expected = K * motion;
    for (Index node = 0; node < 8; ++node) {
        for (Index dof = 0; dof < 6; ++dof) {
            EXPECT_NEAR(rhs(node, dof), expected(6*node+dof), 1e-10);
        }
    }

    // Free compatible expansion must leave no recovered stress or prestress
    model::Field thermal{"THERMAL_FREE_STRAIN", model::FieldDomain::ELEMENT_NODAL,
                         model._data->field_rows(model::FieldDomain::ELEMENT_NODAL), 1};
    thermal.set_zero();
    load.apply_thermal_free_strain(*model._data, thermal);
    const auto recovered = model.compute_stress_nodal(displacement, nullptr, &thermal);
    const auto& stress = std::get<0>(recovered);
    for (Index node = 0; node < 8; ++node) {
        for (Index component = 0; component < 6; ++component) {
            EXPECT_NEAR(stress(node, component), 0.0, 1e-8);
        }
    }
    Precision kg_storage[48*48]{};
    const DynamicMatrix Kg = element->evaluate(nullptr, kg_storage, nullptr, &displacement, nullptr, &thermal, false);
    EXPECT_LT(Kg.norm(), 1e-8);
    model.step_end();
}
