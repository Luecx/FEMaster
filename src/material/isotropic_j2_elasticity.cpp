/**
 * @file isotropic_j2_elasticity.cpp
 * @brief Implements isotropic associative J2 elastoplasticity.
 *
 * The implementation deliberately keeps the public material functions explicit.
 * Persistent state packing, stress reduction, and tangent condensation are written
 * directly where they occur instead of being hidden behind small generic helpers.
 * This makes each constitutive entry point readable in isolation.
 *
 * Internally, the small-strain model uses the classical radial return. The
 * finite-strain model uses
 *
 *     F = Fe Fp
 *
 * with the exponential plastic update
 *
 *     Fp_{n+1} = exp(A) Fp_n.
 *
 * The finite-strain algorithmic tangent is the exact consistent derivative of the
 * converged discrete return map,
 *
 *     C_alg = S_E - S_x R_x^-1 R_E.
 *
 * Every derivative entering this expression is analytical. In particular the
 * matrix exponential is differentiated by its exact spectral Frechet derivative;
 * no constitutive finite differences are used in production code.
 *
 * @author Finn Eggers
 * @date 02.09.2026
 */

#include "isotropic_j2_elasticity.h"

#include "../core/logging.h"
#include "strain/axial_strain_green_lagrange.h"
#include "strain/axial_strain_linearized.h"
#include "strain/shell_material_strain_green_lagrange.h"
#include "strain/shell_material_strain_linearized.h"
#include "strain/volume_strain.h"
#include "strain/volume_strain_green_lagrange.h"
#include "strain/volume_strain_linearized.h"
#include "stress/axial_stress_cauchy.h"
#include "stress/axial_stress_pk2.h"
#include "stress/shell_material_stress_cauchy.h"
#include "stress/shell_material_stress_pk2.h"
#include "stress/volume_stress.h"
#include "stress/volume_stress_cauchy.h"
#include "stress/volume_stress_pk2.h"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace fem::material {
namespace {

constexpr Index state_count = 7;
constexpr Index eqp_index   = 6;

using State      = std::array<Precision, state_count>;
using YieldPoint = IsotropicJ2Elasticity::YieldPoint;
using YieldCurve = std::vector<YieldPoint>;

#include "isotropic_j2_elasticity_utils.ipp"
#include "isotropic_j2_elasticity_small.ipp"
#include "isotropic_j2_elasticity_finite.ipp"
#include "isotropic_j2_elasticity_reductions.ipp"

} // namespace

// -----------------------------------------------------------------------------
// Material definition and capabilities
// -----------------------------------------------------------------------------

/**
 * @brief Constructs an isotropic J2 material from E and nu.
 *
 * The admissible range -1 < nu < 0.5 guarantees positive shear and bulk moduli
 * for positive E.
 */
IsotropicJ2Elasticity::IsotropicJ2Elasticity(Precision youngs_in,
                                             Precision poisson_in)
    : youngs  (youngs_in),
      poisson (poisson_in) {
    logging::error(youngs > Precision(0),
        "J2: Young's modulus must be positive");
    logging::error(poisson > Precision(-1) && poisson < Precision(0.5),
        "J2: Poisson ratio must be in (-1, 0.5)");
}

/**
 * @brief Returns the isotropic shear modulus.
 *
 *     G = E / [2 (1 + nu)].
 */
Precision IsotropicJ2Elasticity::shear_modulus() const {
    return youngs / (Precision(2) * (Precision(1) + poisson));
}

/**
 * @brief Returns the isotropic bulk modulus.
 *
 *     K = E / [3 (1 - 2 nu)].
 */
Precision IsotropicJ2Elasticity::bulk_modulus() const {
    return youngs / (Precision(3) * (Precision(1) - Precision(2) * poisson));
}

/**
 * @brief Appends one point to the tabulated isotropic hardening curve.
 *
 * Equivalent plastic strain must increase strictly. Yield stress must be
 * non-decreasing. The first point is required at zero plastic strain because it
 * defines the initial yield surface.
 */
void IsotropicJ2Elasticity::add_yield_point(Precision yield_stress,
                                            Precision equivalent_plastic_strain) {
    logging::error(std::isfinite(yield_stress) && yield_stress > Precision(0),
        "J2: yield stress must be finite and positive");
    logging::error(std::isfinite(equivalent_plastic_strain)
                   && equivalent_plastic_strain >= Precision(0),
        "J2: equivalent plastic strain must be finite and non-negative");

    const Precision strain_tolerance =
        Precision(100) * std::numeric_limits<Precision>::epsilon();

    if (yield_points_.empty()) {
        logging::error(std::abs(equivalent_plastic_strain) <= strain_tolerance,
            "J2: first yield point must have zero equivalent plastic strain");

        // Remove harmless parser/decimal round-off at the initial point so the
        // stored hardening table starts at exactly alpha = 0.
        equivalent_plastic_strain = Precision(0);
    } else {
        const YieldPoint& previous = yield_points_.back();

        logging::error(equivalent_plastic_strain > previous.equivalent_plastic_strain,
            "J2: equivalent plastic strains must be added in strictly increasing order");
        logging::error(yield_stress >= previous.yield_stress,
            "J2: tabulated isotropic hardening must be non-decreasing");
    }

    yield_points_.push_back({yield_stress, equivalent_plastic_strain});
}

/**
 * @brief Returns the hardening table exactly as supplied to the material.
 */
const std::vector<IsotropicJ2Elasticity::YieldPoint>&
IsotropicJ2Elasticity::get_yield_points() const {
    return yield_points_;
}

/** @brief J2 supports infinitesimal axial strain and Cauchy stress. */
bool IsotropicJ2Elasticity::supports_axial_linearized() const {
    return true;
}

/** @brief J2 supports Green-Lagrange axial strain and PK2 stress. */
bool IsotropicJ2Elasticity::supports_axial_green_lagrange() const {
    return true;
}

/** @brief J2 supports infinitesimal three-dimensional strain. */
bool IsotropicJ2Elasticity::supports_volume_linearized() const {
    return true;
}

/** @brief J2 supports Green-Lagrange three-dimensional strain. */
bool IsotropicJ2Elasticity::supports_volume_green_lagrange() const {
    return true;
}

/** @brief J2 supports infinitesimal integrated-shell material response. */
bool IsotropicJ2Elasticity::supports_shell_integration_linearized() const {
    return true;
}

/** @brief J2 supports finite-strain integrated-shell material response. */
bool IsotropicJ2Elasticity::supports_shell_integration_green_lagrange() const {
    return true;
}

/**
 * @brief Returns the seven persistent J2 history components.
 *
 *     [Cp11, Cp22, Cp33, Cp23, Cp13, Cp12, eqp]
 */
Index IsotropicJ2Elasticity::state_size() const {
    return state_count;
}

/**
 * @brief Initializes the undeformed, virgin material history.
 *
 * Initially
 *
 *     Fp = I
 *     Cp = Fp^T Fp = I
 *     eqp = 0.
 *
 * The explicit assignments intentionally document the persistent state layout at
 * the point where it is initialized.
 */
void IsotropicJ2Elasticity::initialize_state(Precision* state) const {
    state[0] = Precision(1); // Cp11
    state[1] = Precision(1); // Cp22
    state[2] = Precision(1); // Cp33
    state[3] = Precision(0); // Cp23
    state[4] = Precision(0); // Cp13
    state[5] = Precision(0); // Cp12
    state[6] = Precision(0); // accumulated equivalent plastic strain
}

// -----------------------------------------------------------------------------
// Three-dimensional material response
// -----------------------------------------------------------------------------

/**
 * @brief Evaluates the infinitesimal three-dimensional J2 response.
 *
 * `old_state` is copied explicitly into a private working state. The radial return
 * mutates only that local copy. After stress and consistent tangent have been
 * obtained, the complete converged state is written to `new_state` when requested.
 */
void IsotropicJ2Elasticity::evaluate(const VolumeStrainLinearized& strain,
                                     const Precision* old_state,
                                     Precision* new_state,
                                     VolumeStressCauchy& stress,
                                     Mat6& tangent) const {
    // -------------------------------------------------------------------------
    // Start this constitutive candidate from the immutable committed history.
    // -------------------------------------------------------------------------
    State state{
        old_state[0],
        old_state[1],
        old_state[2],
        old_state[3],
        old_state[4],
        old_state[5],
        old_state[6]
    };

    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();

    // VolumeStrainLinearized derives from VolumeStrain, so the internal J2 kernel
    // can use the generic tensor/Voigt representation without knowing the public
    // strain-measure wrapper.
    const SmallResponse response = integrate_small_strain(
        strain,
        state,
        shear,
        bulk,
        yield_points_
    );

    stress  = VolumeStressCauchy(response.stress.voigt());
    tangent = tangent_small(response, shear, bulk, yield_points_);

    // Only the converged candidate history is exposed to the caller. A null
    // new_state is intentionally allowed for state-neutral evaluations.
    if (new_state != nullptr) {
        new_state[0] = state[0];
        new_state[1] = state[1];
        new_state[2] = state[2];
        new_state[3] = state[3];
        new_state[4] = state[4];
        new_state[5] = state[5];
        new_state[6] = state[6];
    }
}

/**
 * @brief Evaluates finite-strain three-dimensional J2 with a required tangent.
 *
 * This overload exists for the normal Elasticity interface and forwards only the
 * optionality of tangent construction to the pointer overload below.
 */
void IsotropicJ2Elasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                     const Precision* old_state,
                                     Precision* new_state,
                                     VolumeStressPK2& stress,
                                     Mat6& tangent) const {
    evaluate(strain, old_state, new_state, stress, &tangent);
}

/**
 * @brief Evaluates the multiplicative finite-strain three-dimensional J2 response.
 *
 * The physical return map is always evaluated. When `tangent == nullptr`, only
 * construction of the expensive consistent tangent is skipped. This is the path
 * used by residual-only global line-search evaluations.
 */
void IsotropicJ2Elasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                     const Precision* old_state,
                                     Precision* new_state,
                                     VolumeStressPK2& stress,
                                     Mat6* tangent) const {
    // Explicit committed -> working state transfer. No previous trial state can
    // enter a new global Newton or line-search candidate.
    State state{
        old_state[0],
        old_state[1],
        old_state[2],
        old_state[3],
        old_state[4],
        old_state[5],
        old_state[6]
    };

    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();

    const FiniteResponse response = integrate_finite_strain(
        strain,
        state,
        shear,
        bulk,
        yield_points_
    );

    stress = VolumeStressPK2(response.stress.voigt());

    if (tangent != nullptr) {
        *tangent = tangent_finite(
            response,
            shear,
            bulk,
            yield_points_
        );
    }

    if (new_state != nullptr) {
        new_state[0] = state[0];
        new_state[1] = state[1];
        new_state[2] = state[2];
        new_state[3] = state[3];
        new_state[4] = state[4];
        new_state[5] = state[5];
        new_state[6] = state[6];
    }
}

// -----------------------------------------------------------------------------
// Shell plane-stress response
// -----------------------------------------------------------------------------

/**
 * @brief Evaluates infinitesimal integrated-shell J2 under sigma_33 = 0.
 *
 * The shell reduction first solves the missing thickness strain with the exact
 * three-dimensional tangent. The returned three-dimensional tangent C is then
 * condensed by the Schur complement
 *
 *     C_shell = C_aa - C_a3 C_3a / C_33,
 *
 * where a denotes the five retained shell strain/stress components.
 */
void IsotropicJ2Elasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                     const Precision* old_state,
                                     Precision* new_state,
                                     ShellMaterialStressCauchy& stress,
                                     Mat5& tangent) const {
    const State committed{
        old_state[0],
        old_state[1],
        old_state[2],
        old_state[3],
        old_state[4],
        old_state[5],
        old_state[6]
    };

    State state{};
    Mat6 tangent_3d;

    const VolumeStress stress_3d = solve_shell_plane_stress<false>(
        strain.values(),
        committed,
        state,
        youngs,
        poisson,
        shear_modulus(),
        bulk_modulus(),
        yield_points_,
        tangent_3d
    );

    // -------------------------------------------------------------------------
    // Map the generic three-dimensional stress back to shell ordering.
    //
    // VolumeStress ordering:
    //     [xx, yy, zz, yz, xz, xy]
    //
    // Shell ordering:
    //     [xx, yy, xy, xz, yz].
    // -------------------------------------------------------------------------
    const Vec6& volume_stress = stress_3d.voigt();
    stress.values() << volume_stress(0),
                       volume_stress(1),
                       volume_stress(5),
                       volume_stress(4),
                       volume_stress(3);

    // -------------------------------------------------------------------------
    // Condense the eliminated zz strain/stress pair from the exact 3D tangent.
    // -------------------------------------------------------------------------
    constexpr std::array<Index, 5> retained{0, 1, 5, 4, 3};
    constexpr Index eliminated = 2;

    Mat5 Caa;
    Eigen::Matrix<Precision, 5, 1> Caz;
    Eigen::Matrix<Precision, 1, 5> Cza;

    for (Index i = 0; i < 5; ++i) {
        const Index row = retained[static_cast<std::size_t>(i)];
        Caz(i) = tangent_3d(row, eliminated);
        Cza(i) = tangent_3d(eliminated, row);

        for (Index j = 0; j < 5; ++j) {
            const Index column = retained[static_cast<std::size_t>(j)];
            Caa(i, j) = tangent_3d(row, column);
        }
    }

    logging::error(std::abs(tangent_3d(eliminated, eliminated))
                   > Precision(100) * std::numeric_limits<Precision>::epsilon(),
        "J2: singular thickness-normal tangent during shell condensation");

    tangent = Caa
            - (Caz * Cza) / tangent_3d(eliminated, eliminated);

    if (new_state != nullptr) {
        new_state[0] = state[0];
        new_state[1] = state[1];
        new_state[2] = state[2];
        new_state[3] = state[3];
        new_state[4] = state[4];
        new_state[5] = state[5];
        new_state[6] = state[6];
    }
}

/**
 * @brief Evaluates finite-strain integrated-shell J2 under S33 = 0.
 *
 * The reduction is identical in structure to the infinitesimal case, but the
 * three-dimensional constitutive kernel uses Green-Lagrange strain, PK2 stress,
 * and the multiplicative finite-strain return map.
 */
void IsotropicJ2Elasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                                     const Precision* old_state,
                                     Precision* new_state,
                                     ShellMaterialStressPK2& stress,
                                     Mat5& tangent) const {
    const State committed{
        old_state[0],
        old_state[1],
        old_state[2],
        old_state[3],
        old_state[4],
        old_state[5],
        old_state[6]
    };

    State state{};
    Mat6 tangent_3d;

    const VolumeStress stress_3d = solve_shell_plane_stress<true>(
        strain.values(),
        committed,
        state,
        youngs,
        poisson,
        shear_modulus(),
        bulk_modulus(),
        yield_points_,
        tangent_3d
    );

    const Vec6& volume_stress = stress_3d.voigt();
    stress.values() << volume_stress(0),
                       volume_stress(1),
                       volume_stress(5),
                       volume_stress(4),
                       volume_stress(3);

    constexpr std::array<Index, 5> retained{0, 1, 5, 4, 3};
    constexpr Index eliminated = 2;

    Mat5 Caa;
    Eigen::Matrix<Precision, 5, 1> Caz;
    Eigen::Matrix<Precision, 1, 5> Cza;

    for (Index i = 0; i < 5; ++i) {
        const Index row = retained[static_cast<std::size_t>(i)];
        Caz(i) = tangent_3d(row, eliminated);
        Cza(i) = tangent_3d(eliminated, row);

        for (Index j = 0; j < 5; ++j) {
            const Index column = retained[static_cast<std::size_t>(j)];
            Caa(i, j) = tangent_3d(row, column);
        }
    }

    logging::error(std::abs(tangent_3d(eliminated, eliminated))
                   > Precision(100) * std::numeric_limits<Precision>::epsilon(),
        "J2: singular thickness-normal tangent during shell condensation");

    tangent = Caa
            - (Caz * Cza) / tangent_3d(eliminated, eliminated);

    if (new_state != nullptr) {
        new_state[0] = state[0];
        new_state[1] = state[1];
        new_state[2] = state[2];
        new_state[3] = state[3];
        new_state[4] = state[4];
        new_state[5] = state[5];
        new_state[6] = state[6];
    }
}

// -----------------------------------------------------------------------------
// Axial uniaxial-stress response
// -----------------------------------------------------------------------------

/**
 * @brief Evaluates infinitesimal axial J2 under sigma_22 = sigma_33 = 0.
 *
 * After solving the two transverse strain components, the scalar tangent is the
 * Schur complement of the transverse 2x2 tangent block,
 *
 *     C_axial = C_11 - C_1b C_bb^-1 C_b1,
 *
 * with b = {22, 33}.
 */
void IsotropicJ2Elasticity::evaluate(const AxialStrainLinearized& strain,
                                     const Precision* old_state,
                                     Precision* new_state,
                                     AxialStressCauchy& stress,
                                     Precision& tangent) const {
    const State committed{
        old_state[0],
        old_state[1],
        old_state[2],
        old_state[3],
        old_state[4],
        old_state[5],
        old_state[6]
    };

    State state{};
    Mat6 tangent_3d;

    const VolumeStress stress_3d = solve_axial_stress<false>(
        strain.value(),
        committed,
        state,
        youngs,
        poisson,
        shear_modulus(),
        bulk_modulus(),
        yield_points_,
        tangent_3d
    );

    stress.value() = stress_3d[VolumeStress::Component::XX];

    Eigen::Matrix<Precision, 1, 2> C1b;
    Eigen::Matrix<Precision, 2, 1> Cb1;
    Eigen::Matrix<Precision, 2, 2> Cbb;

    C1b << tangent_3d(0, 1), tangent_3d(0, 2);
    Cb1 << tangent_3d(1, 0), tangent_3d(2, 0);
    Cbb << tangent_3d(1, 1), tangent_3d(1, 2),
           tangent_3d(2, 1), tangent_3d(2, 2);

    Eigen::FullPivLU<Eigen::Matrix<Precision, 2, 2>> lu(Cbb);
    logging::error(lu.isInvertible(),
        "J2: singular transverse block during axial tangent condensation");

    tangent = tangent_3d(0, 0) - (C1b * lu.solve(Cb1))(0, 0);

    if (new_state != nullptr) {
        new_state[0] = state[0];
        new_state[1] = state[1];
        new_state[2] = state[2];
        new_state[3] = state[3];
        new_state[4] = state[4];
        new_state[5] = state[5];
        new_state[6] = state[6];
    }
}

/**
 * @brief Evaluates finite-strain axial J2 under S22 = S33 = 0.
 *
 * Green-Lagrange axial strain and PK2 stress are work-conjugate. The transverse
 * stress solve and final scalar Schur complement therefore use the exact finite-
 * strain PK2/Green-Lagrange tangent without any change of stress measure.
 */
void IsotropicJ2Elasticity::evaluate(const AxialStrainGreenLagrange& strain,
                                     const Precision* old_state,
                                     Precision* new_state,
                                     AxialStressPK2& stress,
                                     Precision& tangent) const {
    const State committed{
        old_state[0],
        old_state[1],
        old_state[2],
        old_state[3],
        old_state[4],
        old_state[5],
        old_state[6]
    };

    State state{};
    Mat6 tangent_3d;

    const VolumeStress stress_3d = solve_axial_stress<true>(
        strain.value(),
        committed,
        state,
        youngs,
        poisson,
        shear_modulus(),
        bulk_modulus(),
        yield_points_,
        tangent_3d
    );

    stress.value() = stress_3d[VolumeStress::Component::XX];

    Eigen::Matrix<Precision, 1, 2> C1b;
    Eigen::Matrix<Precision, 2, 1> Cb1;
    Eigen::Matrix<Precision, 2, 2> Cbb;

    C1b << tangent_3d(0, 1), tangent_3d(0, 2);
    Cb1 << tangent_3d(1, 0), tangent_3d(2, 0);
    Cbb << tangent_3d(1, 1), tangent_3d(1, 2),
           tangent_3d(2, 1), tangent_3d(2, 2);

    Eigen::FullPivLU<Eigen::Matrix<Precision, 2, 2>> lu(Cbb);
    logging::error(lu.isInvertible(),
        "J2: singular transverse block during axial tangent condensation");

    tangent = tangent_3d(0, 0) - (C1b * lu.solve(Cb1))(0, 0);

    if (new_state != nullptr) {
        new_state[0] = state[0];
        new_state[1] = state[1];
        new_state[2] = state[2];
        new_state[3] = state[3];
        new_state[4] = state[4];
        new_state[5] = state[5];
        new_state[6] = state[6];
    }
}

} // namespace fem::material
