/**
 * @file isotropic_j2_elasticity.cpp
 * @brief Implements isotropic associative J2 elastoplasticity.
 *
 * The implementation deliberately keeps the public material functions explicit.
 * Persistent state packing, stress reduction and tangent condensation are written
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
 *     Fp_(n+1) = exp(A) Fp_n.
 *
 * The finite-strain algorithmic tangent is the exact consistent derivative of
 * the converged discrete return map,
 *
 *     C_alg = S_E - S_x R_x^-1 R_E.
 *
 * Every derivative entering this expression is analytical. In particular the
 * matrix exponential is differentiated by its exact spectral Frechet derivative;
 * no constitutive finite differences are used in production code.
 *
 * Tangent output follows the nullable-pointer contract of `Elasticity`. The full
 * three-dimensional constitutive paths therefore omit tangent construction when
 * the caller requests stress and state only. Axial and shell reductions still
 * form a three-dimensional tangent internally because their local zero-stress
 * Newton equations require that derivative.
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
 * Constructs isotropic J2 plasticity from Young's modulus and Poisson's ratio.
 *
 * Positive `E` and the open interval
 *
 *     -1 < nu < 0.5
 *
 * guarantee positive shear and bulk moduli. Plastic hardening data is appended
 * separately with `add_yield_point()`.
 *
 * @param youngs_in Young's modulus E.
 * @param poisson_in Poisson ratio nu.
 */
IsotropicJ2Elasticity::IsotropicJ2Elasticity(Precision youngs_in,
                                             Precision poisson_in)
    : youngs (youngs_in),
      poisson(poisson_in) {
    logging::error(youngs > Precision(0),
        "J2: Young's modulus must be positive");
    logging::error(poisson > Precision(-1) && poisson < Precision(0.5),
        "J2: Poisson ratio must be in (-1, 0.5)");
}

/**
 * Returns the elastic shear modulus.
 *
 *     G = E / [2 (1 + nu)].
 *
 * @return Elastic shear modulus G.
 */
Precision IsotropicJ2Elasticity::shear_modulus() const {
    return youngs / (Precision(2) * (Precision(1) + poisson));
}

/**
 * Returns the elastic bulk modulus.
 *
 *     K = E / [3 (1 - 2 nu)].
 *
 * @return Elastic bulk modulus K.
 */
Precision IsotropicJ2Elasticity::bulk_modulus() const {
    return youngs / (Precision(3) * (Precision(1) - Precision(2) * poisson));
}

/**
 * Appends one point to the tabulated isotropic hardening curve.
 *
 * Equivalent plastic strain must increase strictly and yield stress must be
 * non-decreasing. The first point is required at zero plastic strain because it
 * defines the initial yield surface.
 *
 * @param yield_stress Yield stress at the supplied plastic strain.
 * @param equivalent_plastic_strain Accumulated equivalent plastic strain.
 */
void IsotropicJ2Elasticity::add_yield_point(Precision yield_stress,
                                            Precision equivalent_plastic_strain) {
    // Validate the physical range before touching the hardening table.
    logging::error(std::isfinite(yield_stress) && yield_stress > Precision(0),
        "J2: yield stress must be finite and positive");
    logging::error(std::isfinite(equivalent_plastic_strain)
                   && equivalent_plastic_strain >= Precision(0),
        "J2: equivalent plastic strain must be finite and non-negative");

    const Precision strain_tolerance =
        Precision(100) * std::numeric_limits<Precision>::epsilon();

    if (yield_points_.empty()) {
        // The first point defines initial yield and therefore belongs at alpha = 0.
        logging::error(std::abs(equivalent_plastic_strain) <= strain_tolerance,
            "J2: first yield point must have zero equivalent plastic strain");

        // Remove harmless parser/decimal round-off from that exact initial value.
        equivalent_plastic_strain = Precision(0);
    } else {
        // Every later point must extend the existing piecewise-linear curve.
        const YieldPoint& previous = yield_points_.back();

        logging::error(equivalent_plastic_strain > previous.equivalent_plastic_strain,
            "J2: equivalent plastic strains must be added in strictly increasing order");
        logging::error(yield_stress >= previous.yield_stress,
            "J2: tabulated isotropic hardening must be non-decreasing");
    }

    yield_points_.push_back({yield_stress, equivalent_plastic_strain});
}

/**
 * Returns the hardening table exactly as supplied to the material.
 *
 * @return Piecewise-linear isotropic hardening points in ascending plastic strain.
 */
const std::vector<IsotropicJ2Elasticity::YieldPoint>&
IsotropicJ2Elasticity::get_yield_points() const {
    return yield_points_;
}

bool IsotropicJ2Elasticity::supports_axial_linearized() const {
    return true;
}

bool IsotropicJ2Elasticity::supports_axial_green_lagrange() const {
    return true;
}

bool IsotropicJ2Elasticity::supports_volume_linearized() const {
    return true;
}

bool IsotropicJ2Elasticity::supports_volume_green_lagrange() const {
    return true;
}

bool IsotropicJ2Elasticity::supports_shell_integration_linearized() const {
    return true;
}

bool IsotropicJ2Elasticity::supports_shell_integration_green_lagrange() const {
    return true;
}

/**
 * Returns the number of persistent J2 history components.
 *
 * The state layout is
 *
 *     [Cp11, Cp22, Cp33, Cp23, Cp13, Cp12, eqp].
 *
 * @return Seven persistent scalar state components.
 */
Index IsotropicJ2Elasticity::state_size() const {
    return state_count;
}

/**
 * Initializes the undeformed, virgin material history.
 *
 * Initially
 *
 *     Fp  = I
 *     Cp  = Fp^T Fp
 *         = I
 *     eqp = 0.
 *
 * The explicit assignments document the persistent state layout at the point
 * where the state is created.
 *
 * @param state Writable material-point state row.
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
 * Evaluates the infinitesimal three-dimensional J2 response.
 *
 * Every constitutive candidate starts from the immutable committed state. The
 * small-strain radial return mutates only a local working copy. Cauchy stress is
 * always returned, while the exact radial-return tangent is constructed only
 * when the caller supplies tangent storage.
 *
 * After the constitutive update, the converged candidate history is copied to
 * `new_state` when requested. Passing `new_state == nullptr` therefore performs
 * an otherwise identical state-neutral material evaluation.
 *
 * @param strain Infinitesimal volume strain in the material basis.
 * @param old_state Immutable committed J2 state.
 * @param new_state Optional destination for the converged trial state.
 * @param stress Cauchy stress in the material basis.
 * @param tangent Optional consistent derivative `d sigma/d epsilon`.
 */
void IsotropicJ2Elasticity::evaluate(const VolumeStrainLinearized& strain,
                                     const Precision*              old_state,
                                     Precision*                    new_state,
                                     VolumeStressCauchy&           stress,
                                     Mat6*                         tangent) const {
    // Start the candidate from the committed history. Writing all seven entries
    // explicitly keeps the state layout visible and prevents any previous trial
    // state from feeding the next Newton or line-search evaluation.
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

    // VolumeStrainLinearized derives from VolumeStrain, so the internal return
    // map works entirely with the common six-component volume representation.
    const SmallResponse response = integrate_small_strain(
        strain,
        state,
        shear,
        bulk,
        yield_points_
    );

    // Stress belongs to every constitutive evaluation. The tangent is additional
    // output and is intentionally omitted for stress-only calls.
    stress = VolumeStressCauchy(response.stress.voigt());

    if (tangent != nullptr) {
        *tangent = tangent_small(response, shear, bulk, yield_points_);
    }

    // Publish the converged candidate state only when persistent trial storage was
    // supplied by the nonlinear state manager.
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
 * Evaluates the multiplicative finite-strain three-dimensional J2 response.
 *
 * Green-Lagrange strain supplies the total metric
 *
 *     C = I + 2 E.
 *
 * The finite return map integrates the plastic metric and returns PK2 stress.
 * The algorithmic tangent is the consistent analytic derivative of exactly that
 * converged return map and is evaluated only for a non-null tangent pointer.
 * This is the relevant fast path for residual-only global line-search assembly.
 *
 * @param strain Green-Lagrange strain in the reference material basis.
 * @param old_state Immutable committed J2 state.
 * @param new_state Optional destination for the converged trial state.
 * @param stress Second Piola-Kirchhoff stress.
 * @param tangent Optional consistent derivative `dS/dE`.
 */
void IsotropicJ2Elasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                     const Precision*                 old_state,
                                     Precision*                       new_state,
                                     VolumeStressPK2&                 stress,
                                     Mat6*                            tangent) const {
    // Reconstruct a private working state directly from the committed row.
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

    // Integrate the complete physical state regardless of whether the caller also
    // requests a tangent. This preserves identical stress/state semantics between
    // Newton and residual-only line-search evaluations.
    const FiniteResponse response = integrate_finite_strain(
        strain,
        state,
        shear,
        bulk,
        yield_points_
    );

    stress = VolumeStressPK2(response.stress.voigt());

    // The consistent finite-strain tangent is the expensive part of the response,
    // so construct it only when global assembly actually needs it.
    if (tangent != nullptr) {
        *tangent = tangent_finite(response, shear, bulk, yield_points_);
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
 * Evaluates infinitesimal integrated-shell J2 under `sigma33 = 0`.
 *
 * The shell strain contains five retained components. The missing thickness
 * strain is determined by a local Newton solve using the exact three-dimensional
 * J2 tangent. That local derivative is required even when the external shell
 * tangent is not requested.
 *
 * When tangent output is requested, the converged three-dimensional tangent is
 * condensed by
 *
 *     C_shell = C_aa - C_a3 C_3a / C_33,
 *
 * where `a = {11,22,12,13,23}` denotes the retained shell components.
 *
 * @param strain Linearized shell material strain.
 * @param old_state Immutable committed J2 state.
 * @param new_state Optional destination for the converged trial state.
 * @param stress Plane-stress shell Cauchy stress.
 * @param tangent Optional condensed shell tangent.
 */
void IsotropicJ2Elasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                     const Precision*                     old_state,
                                     Precision*                           new_state,
                                     ShellMaterialStressCauchy&            stress,
                                     Mat5*                                tangent) const {
    // Preserve the committed state while the local thickness-strain Newton solve
    // evaluates several constitutive candidates.
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
    Mat6  tangent_3d;

    // Solve sigma33 = 0. The reduction itself needs tangent_3d for its Newton
    // derivative independently of whether the caller wants the final Mat5.
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

    // Map volume ordering [xx,yy,zz,yz,xz,xy] to shell ordering
    // [xx,yy,xy,xz,yz].
    const Vec6& volume_stress = stress_3d.voigt();
    stress.values() << volume_stress(0),
                       volume_stress(1),
                       volume_stress(5),
                       volume_stress(4),
                       volume_stress(3);

    if (tangent != nullptr) {
        // Extract the retained blocks around the eliminated zz component.
        constexpr std::array<Index, 5> retained {0, 1, 5, 4, 3};
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

        *tangent = Caa - (Caz * Cza) / tangent_3d(eliminated, eliminated);
    }

    // Publish only the state belonging to the converged plane-stress candidate.
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
 * Evaluates finite-strain integrated-shell J2 under `S33 = 0`.
 *
 * The reduction has the same structure as the infinitesimal shell case, but the
 * three-dimensional constitutive candidates use Green-Lagrange strain, PK2
 * stress and the multiplicative finite-strain return map. The local thickness
 * solve always requires the three-dimensional consistent tangent; the final
 * five-component Schur complement is produced only when requested.
 *
 * @param strain Green-Lagrange shell material strain.
 * @param old_state Immutable committed J2 state.
 * @param new_state Optional destination for the converged trial state.
 * @param stress Plane-stress shell PK2 stress.
 * @param tangent Optional condensed consistent shell tangent.
 */
void IsotropicJ2Elasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                                     const Precision*                        old_state,
                                     Precision*                              new_state,
                                     ShellMaterialStressPK2&                 stress,
                                     Mat5*                                   tangent) const {
    // Keep every local thickness-strain candidate anchored to the same committed
    // history. Only the converged candidate is exposed at the end.
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
    Mat6  tangent_3d;

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

    // Convert the generic volume stress into the shell material component order.
    const Vec6& volume_stress = stress_3d.voigt();
    stress.values() << volume_stress(0),
                       volume_stress(1),
                       volume_stress(5),
                       volume_stress(4),
                       volume_stress(3);

    if (tangent != nullptr) {
        // Condense the eliminated thickness-normal pair from the converged exact
        // three-dimensional PK2/Green-Lagrange tangent.
        constexpr std::array<Index, 5> retained {0, 1, 5, 4, 3};
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

        *tangent = Caa - (Caz * Cza) / tangent_3d(eliminated, eliminated);
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
// Axial uniaxial-stress response
// -----------------------------------------------------------------------------

/**
 * Evaluates infinitesimal axial J2 under `sigma22 = sigma33 = 0`.
 *
 * The two missing transverse strains are determined by the three-dimensional
 * constitutive reduction. Its local Newton solve always requires the full
 * tangent. When the external scalar tangent is requested, the transverse block
 * is eliminated by
 *
 *     C_axial = C_11 - C_1b C_bb^-1 C_b1,
 *
 * with `b = {22,33}`.
 *
 * @param strain Infinitesimal axial strain.
 * @param old_state Immutable committed J2 state.
 * @param new_state Optional destination for the converged trial state.
 * @param stress Axial Cauchy stress.
 * @param tangent Optional consistent scalar axial tangent.
 */
void IsotropicJ2Elasticity::evaluate(const AxialStrainLinearized& strain,
                                     const Precision*             old_state,
                                     Precision*                   new_state,
                                     AxialStressCauchy&           stress,
                                     Precision*                   tangent) const {
    // Keep the transverse local solve anchored to the committed material history.
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
    Mat6  tangent_3d;

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

    // Return the retained axial stress component from the converged 3D state.
    stress.value() = stress_3d[VolumeStress::Component::XX];

    if (tangent != nullptr) {
        // Extract the transverse coupling blocks needed by the Schur complement.
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

        *tangent = tangent_3d(0, 0) - (C1b * lu.solve(Cb1))(0, 0);
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

/**
 * Evaluates finite-strain axial J2 under `S22 = S33 = 0`.
 *
 * Green-Lagrange axial strain and PK2 stress are work-conjugate. The transverse
 * stress solve and optional scalar Schur complement therefore use the exact
 * finite-strain PK2/Green-Lagrange tangent without changing stress measure.
 * The three-dimensional tangent remains necessary internally for the transverse
 * Newton solve even when the final scalar tangent is not requested.
 *
 * @param strain Axial Green-Lagrange strain.
 * @param old_state Immutable committed J2 state.
 * @param new_state Optional destination for the converged trial state.
 * @param stress Axial second Piola-Kirchhoff stress.
 * @param tangent Optional consistent scalar derivative `dS11/dE11`.
 */
void IsotropicJ2Elasticity::evaluate(const AxialStrainGreenLagrange& strain,
                                     const Precision*                old_state,
                                     Precision*                      new_state,
                                     AxialStressPK2&                 stress,
                                     Precision*                      tangent) const {
    // Every transverse Newton candidate starts from the same committed state.
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
    Mat6  tangent_3d;

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

    if (tangent != nullptr) {
        // Condense the two traction-free transverse components from the exact
        // three-dimensional algorithmic tangent.
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

        *tangent = tangent_3d(0, 0) - (C1b * lu.solve(Cb1))(0, 0);
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

} // namespace fem::material
