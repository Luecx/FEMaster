/**
 * @file test_material_j2.cpp
 * @brief Verifies J2 stress updates, plastic history and material tangents.
 *
 * The tests exercise small-strain and multiplicative finite-strain return maps,
 * axial reductions and parser registration. Near-incompressible cases verify
 * the physical Cauchy yield surface and volume-preserving plastic history.
 */

#include "../src/io/dsl/deck_parser.h"
#include "../src/io/dsl/file.h"
#include "../src/io/dsl/registry.h"
#include "../src/io/reader/commands/register_functions.h"
#include "../src/material/isotropic_j2_elasticity.h"
#include "../src/material/strain/axial_strain_green_lagrange.h"
#include "../src/material/strain/volume_strain_green_lagrange.h"
#include "../src/material/strain/volume_strain_linearized.h"
#include "../src/material/stress/volume_stress_cauchy.h"
#include "../src/material/stress/axial_stress_pk2.h"
#include "../src/material/stress/volume_stress_pk2.h"
#include "../src/model/model.h"

#include <gtest/gtest.h>
#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace fem;

TEST(Material_J2, DerivedModuliAndYieldPoints) {
    material::IsotropicJ2Elasticity j2(Precision(210000), Precision(0.3));

    EXPECT_NEAR(j2.shear_modulus(), Precision(80769.23076923077), Precision(1e-8));
    EXPECT_NEAR(j2.bulk_modulus(), Precision(175000), Precision(1e-8));
    EXPECT_EQ(j2.state_size(), 7);

    j2.add_yield_point(Precision(280), Precision(0));
    j2.add_yield_point(Precision(320), Precision(0.02));
    j2.add_yield_point(Precision(380), Precision(0.10));

    const auto& points = j2.get_yield_points();
    ASSERT_EQ(points.size(), 3u);
    EXPECT_EQ(points[0].yield_stress, Precision(280));
    EXPECT_EQ(points[0].equivalent_plastic_strain, Precision(0));
    EXPECT_EQ(points[1].yield_stress, Precision(320));
    EXPECT_EQ(points[1].equivalent_plastic_strain, Precision(0.02));
    EXPECT_EQ(points[2].yield_stress, Precision(380));
    EXPECT_EQ(points[2].equivalent_plastic_strain, Precision(0.10));
}

TEST(Material_J2, EvaluationAlwaysStartsFromCommittedState) {
    material::IsotropicJ2Elasticity j2(Precision(210000), Precision(0.3));
    j2.add_yield_point(Precision(280), Precision(0));
    j2.add_yield_point(Precision(350), Precision(0.10));

    std::vector<Precision> committed(static_cast<std::size_t>(j2.state_size()));
    j2.initialize_state(committed.data());
    const auto committed_before = committed;

    // A sufficiently large deviatoric strain drives the material into plasticity.
    Vec6 strain_values = Vec6::Zero();
    strain_values(0) = Precision(0.01);
    const VolumeStrainLinearized strain(strain_values);

    std::vector<Precision> trial_a(committed.size());
    std::vector<Precision> trial_b(committed.size());
    VolumeStressCauchy stress_a;
    VolumeStressCauchy stress_b;
    VolumeStressCauchy stress_readonly;
    Mat6 tangent_a;
    Mat6 tangent_b;
    Mat6 tangent_readonly;

    j2.evaluate(strain, committed.data(), trial_a.data(), stress_a, &tangent_a);
    j2.evaluate(strain, committed.data(), trial_b.data(), stress_b, &tangent_b);
    j2.evaluate(strain, committed.data(), nullptr, stress_readonly, &tangent_readonly);

    // Constitutive calls must never modify committed history and repeated trial
    // evaluations from the same committed state must be deterministic.
    for (std::size_t i = 0; i < committed.size(); ++i) {
        EXPECT_EQ(committed[i], committed_before[i]);
        EXPECT_NEAR(trial_a[i], trial_b[i], Precision(1e-12));
    }

    EXPECT_GT(trial_a[6], Precision(0));
    EXPECT_TRUE(stress_a.voigt().isApprox(stress_b.voigt(), Precision(1e-10)));
    EXPECT_TRUE(stress_a.voigt().isApprox(stress_readonly.voigt(), Precision(1e-10)));
    EXPECT_TRUE(tangent_a.isApprox(tangent_b, Precision(1e-10)));
    EXPECT_TRUE(tangent_a.isApprox(tangent_readonly, Precision(1e-10)));
}

TEST(Material_J2, NearIncompressibleYieldCheckUsesDeviatoricScale) {
    // This Poisson ratio is intentionally extremely close to the incompressible
    // limit. The bulk modulus is therefore enormous while the shear modulus and
    // J2 yield surface remain perfectly regular. A yield tolerance scaled with K
    // would incorrectly classify the trial state below as elastic.
    material::IsotropicJ2Elasticity j2(
        Precision(210000),
        Precision(0.499999999)
    );
    j2.add_yield_point(Precision(250), Precision(0));
    j2.add_yield_point(Precision(500), Precision(0.10));

    // Use a trace-free normal strain so the trial response is purely deviatoric.
    // For G ~= 70000 MPa and e = diag(a,-a,0), a = 0.005 produces a trial
    // equivalent stress well above the 250 MPa yield stress, independently of K.
    Vec6 strain_values = Vec6::Zero();
    strain_values(0) = Precision(0.005);
    strain_values(1) = Precision(-0.005);

    std::vector<Precision> committed(static_cast<std::size_t>(j2.state_size()));
    std::vector<Precision> trial_small(committed.size());
    std::vector<Precision> trial_finite(committed.size());
    j2.initialize_state(committed.data());

    VolumeStressCauchy small_stress;
    j2.evaluate(
        VolumeStrainLinearized(strain_values),
        committed.data(),
        trial_small.data(),
        small_stress,
        nullptr
    );

    VolumeStressPK2 finite_stress;
    j2.evaluate(
        VolumeStrainGreenLagrange(strain_values),
        committed.data(),
        trial_finite.data(),
        finite_stress,
        nullptr
    );

    // Both formulations must enter plasticity. In the J2 state layout the final
    // component is accumulated equivalent plastic strain, so a positive value is
    // an unambiguous regression check for the yield decision.
    EXPECT_GT(trial_small[6], Precision(0));
    EXPECT_GT(trial_finite[6], Precision(0));
}

TEST(Material_J2, NearIncompressibleFiniteReturnSatisfiesRotatedYieldSurface) {
    // Rotate a trace-free Green-Lagrange strain to exercise off-diagonal
    // metrics as well as the principal-axis case at several bulk/shear ratios.
    for (const Precision poisson : {Precision(0.49), Precision(0.4999), Precision(0.499999), Precision(0.499999999)}) {
        for (const Precision angle : {Precision(0), Precision(0.7)}) {
            SCOPED_TRACE(::testing::Message() << "nu=" << poisson << " angle=" << angle);
            material::IsotropicJ2Elasticity j2(Precision(210000), poisson);
            j2.add_yield_point(Precision(250), Precision(0));
            j2.add_yield_point(Precision(500), Precision(0.1));
            std::vector<Precision> committed(j2.state_size());
            std::vector<Precision> trial(j2.state_size());
            j2.initialize_state(committed.data());

            const Mat3 rotation = Eigen::AngleAxis<Precision>(angle, Vec3(1, 2, 3).normalized()).toRotationMatrix();
            Mat3 principal = Mat3::Zero();
            principal(0, 0) = Precision(0.005);
            principal(1, 1) = Precision(-0.005);
            const Mat3 strain = rotation * principal * rotation.transpose();
            VolumeStressPK2 stress;
            Mat6 tangent;
            j2.evaluate(VolumeStrainGreenLagrange(strain), committed.data(), trial.data(), stress, &tangent);
            ASSERT_GT(trial[6], Precision(0));
            EXPECT_TRUE(tangent.allFinite());
            EXPECT_TRUE(tangent.isApprox(tangent.transpose(), Precision(1e-12)));

            // Push PK2 forward with F=sqrt(C) and check the tabulated Cauchy
            // yield radius independently of the internal Mandel-stress residual.
            Mat3 stretch = Mat3::Identity();
            stretch(0, 0) = std::sqrt(Precision(1.01));
            stretch(1, 1) = std::sqrt(Precision(0.99));
            const Mat3 F      = rotation * stretch * rotation.transpose();
            const Mat3 cauchy = stress.to_cauchy(F).tensor();
            const Mat3 dev    = cauchy - cauchy.trace() / Precision(3) * Mat3::Identity();
            const Precision q = std::sqrt(Precision(1.5) * dev.squaredNorm());
            EXPECT_NEAR(q, Precision(250) + Precision(2500) * trial[6], Precision(1e-4));
            Mat3 Cp;
            Cp << trial[0], trial[5], trial[4],
                  trial[5], trial[1], trial[3],
                  trial[4], trial[3], trial[2];
            EXPECT_NEAR(Cp.determinant(), Precision(1), Precision(1e-12));
        }
    }
}

TEST(Material_J2, FiniteStrainPlasticTableUsesCauchyYieldStress) {
    material::IsotropicJ2Elasticity j2(Precision(210000), Precision(0.3));
    j2.add_yield_point(Precision(250), Precision(0));

    std::vector<Precision> committed(static_cast<std::size_t>(j2.state_size()));
    std::vector<Precision> trial(committed.size());
    j2.initialize_state(committed.data());

    // Add a moderate volumetric expansion and a small deviatoric perturbation.
    // For the virgin finite-strain trial state this gives approximately
    //
    //     J       = 1.15369
    //     q_M     = 269.47 MPa
    //     J sigma = 288.42 MPa.
    //
    // The raw table value 250 MPa would incorrectly trigger yielding if it were
    // compared directly with Mandel stress. As a true/Cauchy yield stress it must
    // first be converted to the Kirchhoff/Mandel radius with J.
    Vec6 strain_values = Vec6::Zero();
    strain_values(0) = Precision(0.05078);
    strain_values(1) = Precision(0.04961);
    strain_values(2) = Precision(0.04961);

    VolumeStressPK2 stress;
    j2.evaluate(
        VolumeStrainGreenLagrange(strain_values),
        committed.data(),
        trial.data(),
        stress,
        nullptr
    );

    EXPECT_NEAR(trial[6], Precision(0), Precision(1e-14));
}

TEST(Material_J2, FiniteStrainTangentIsSymmetricApproximationOfReturnMap) {
    material::IsotropicJ2Elasticity j2(Precision(210000), Precision(0.3));

    // A single linear hardening segment avoids derivative ambiguity at interior
    // table knots while still exercising a genuinely plastic committed state.
    j2.add_yield_point(Precision(250), Precision(0));
    j2.add_yield_point(Precision(700), Precision(0.30));

    std::vector<Precision> initial(static_cast<std::size_t>(j2.state_size()));
    std::vector<Precision> committed(initial.size());
    j2.initialize_state(initial.data());

    // First create non-trivial committed plastic history.
    Vec6 preload_values = Vec6::Zero();
    preload_values << Precision(0.012), Precision(-0.002), Precision(-0.001),
                      Precision(0.0005), Precision(-0.0003), Precision(0.0004);

    VolumeStressPK2 preload_stress;
    Mat6 preload_tangent;
    j2.evaluate(
        VolumeStrainGreenLagrange(preload_values),
        initial.data(),
        committed.data(),
        preload_stress,
        &preload_tangent
    );
    ASSERT_GT(committed[6], Precision(0));

    // Evaluate the solver tangent at a nearby plastic state from exactly the
    // same committed history used by every independent perturbation below.
    Vec6 strain_values = preload_values;
    strain_values(0) += Precision(0.003);
    strain_values(1) -= Precision(0.0004);
    strain_values(3) += Precision(0.0002);

    VolumeStressPK2 stress;
    Mat6 tangent;
    j2.evaluate(
        VolumeStrainGreenLagrange(strain_values),
        committed.data(),
        nullptr,
        stress,
        &tangent
    );

    // Finite differences recover the exact derivative of the unchanged stress
    // update. The solver tangent is deliberately its symmetric part.
    Mat6 tangent_fd = Mat6::Zero();
    const Precision h = Precision(1e-7);

    for (Index column = 0; column < 6; ++column) {
        Vec6 plus  = strain_values;
        Vec6 minus = strain_values;
        plus(column)  += h;
        minus(column) -= h;

        VolumeStressPK2 stress_plus;
        VolumeStressPK2 stress_minus;
        j2.evaluate(
            VolumeStrainGreenLagrange(plus),
            committed.data(),
            nullptr,
            stress_plus,
            nullptr
        );
        j2.evaluate(
            VolumeStrainGreenLagrange(minus),
            committed.data(),
            nullptr,
            stress_minus,
            nullptr
        );

        tangent_fd.col(column) =
            (stress_plus.voigt() - stress_minus.voigt()) / (Precision(2) * h);
    }

    const Mat6 tangent_fd_symmetric =
        Precision(0.5) * (tangent_fd + tangent_fd.transpose());
    const Precision scale = std::max(Precision(1), tangent_fd.norm());

    EXPECT_LT((tangent - tangent.transpose()).norm() / scale, Precision(1e-12));
    EXPECT_LT((tangent - tangent_fd_symmetric).norm() / scale, Precision(2e-4));

    // PLASTIC remains a true/Cauchy hardening table and the stress/state update
    // remains unchanged. Only the Newton tangent is symmetrized so the global
    // nonlinear solve can use a symmetric-indefinite sparse factorization.
}

TEST(Material_J2, FiniteAxialTangentRemainsConsistentWithReducedStress) {
    material::IsotropicJ2Elasticity j2(Precision(210000), Precision(0.3));
    j2.add_yield_point(Precision(250), Precision(0));
    j2.add_yield_point(Precision(700), Precision(0.30));

    std::vector<Precision> initial(static_cast<std::size_t>(j2.state_size()));
    std::vector<Precision> committed(initial.size());
    j2.initialize_state(initial.data());

    AxialStressPK2 preload_stress;
    Precision preload_tangent = Precision(0);
    j2.evaluate(
        AxialStrainGreenLagrange(Precision(0.012)),
        initial.data(),
        committed.data(),
        preload_stress,
        &preload_tangent
    );
    ASSERT_GT(committed[6], Precision(0));

    const Precision strain_value = Precision(0.015);
    AxialStressPK2 stress;
    Precision tangent = Precision(0);
    j2.evaluate(
        AxialStrainGreenLagrange(strain_value),
        committed.data(),
        nullptr,
        stress,
        &tangent
    );

    const Precision h = Precision(1e-7);
    AxialStressPK2 stress_plus;
    AxialStressPK2 stress_minus;

    j2.evaluate(
        AxialStrainGreenLagrange(strain_value + h),
        committed.data(),
        nullptr,
        stress_plus,
        nullptr
    );
    j2.evaluate(
        AxialStrainGreenLagrange(strain_value - h),
        committed.data(),
        nullptr,
        stress_minus,
        nullptr
    );

    const Precision tangent_fd =
        (stress_plus.value() - stress_minus.value()) / (Precision(2) * h);
    const Precision scale = std::max(Precision(1), std::abs(tangent_fd));

    EXPECT_LT(std::abs(tangent - tangent_fd) / scale, Precision(2e-4));
}

TEST(Material_J2, PlasticCommandReplacesIsotropicElasticity) {
    const std::string input_path = "tests/TMP_J2_MATERIAL.INP";
    std::filesystem::remove(input_path);

    // Intentionally place PLASTIC before ELASTIC in the source. The semantic
    // material order must still execute ELASTIC first and PLASTIC afterwards.
    {
        std::ofstream os(input_path);
        ASSERT_TRUE(os.is_open());
        os << "*MATERIAL, NAME=STEEL\n";
        os << "*PLASTIC\n";
        os << "280., 0.\n";
        os << "320., 0.02\n";
        os << "380., 0.10\n";
        os << "*ELASTIC\n";
        os << "210000., 0.3\n";
    }

    model::Model model;
    io::dsl::Registry registry;
    io::reader::commands::register_material(registry, model);
    io::reader::commands::register_elastic(registry, model);
    io::reader::commands::register_plastic(registry, model);

    io::dsl::File file(input_path);
    io::dsl::DeckParser parser(registry);
    const io::dsl::Deck deck = parser.parse(file);

    const auto materials = deck.root().children("MATERIAL");
    ASSERT_EQ(materials.size(), 1u);

    const auto* definition = materials.front();
    definition->enter();
    definition->execute_children("ELASTIC");
    definition->execute_children("PLASTIC");
    definition->leave();

    const auto steel = model._data->materials.get("STEEL");
    ASSERT_NE(steel, nullptr);
    ASSERT_TRUE(steel->has_elasticity());

    const auto* j2 = steel->elasticity()->as<material::IsotropicJ2Elasticity>();
    ASSERT_NE(j2, nullptr);
    EXPECT_EQ(j2->youngs, Precision(210000));
    EXPECT_EQ(j2->poisson, Precision(0.3));

    const auto& points = j2->get_yield_points();
    ASSERT_EQ(points.size(), 3u);
    EXPECT_EQ(points[0].yield_stress, Precision(280));
    EXPECT_EQ(points[1].yield_stress, Precision(320));
    EXPECT_EQ(points[2].yield_stress, Precision(380));

    std::filesystem::remove(input_path);
}
