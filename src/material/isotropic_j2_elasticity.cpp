/**
 * @file isotropic_j2_elasticity.cpp
 * @brief Implements isotropic small- and finite-strain J2 elastoplasticity.
 *
 * Finite strain follows the multiplicative split F = Fe Fp without requiring F
 * in the public constitutive interface. For an isotropic material the plastic
 * rotation is constitutively irrelevant, so the history is represented by the
 * plastic metric Cp = Fp^T Fp. The total metric follows exactly from the supplied
 * Green-Lagrange strain: C = I + 2 E.
 *
 * The elastic part intentionally matches IsotropicElasticity. In the intermediate
 * configuration
 *
 *     Ee = 1/2 (Ce - I),
 *     Se = lambda tr(Ee) I + 2 G Ee,
 *     M  = Ce Se.
 *
 * Thus an initially elastic finite-strain response is exactly the same
 * Saint-Venant-Kirchhoff law used by IsotropicElasticity. Plastic yielding uses
 *
 *     q = sqrt(3/2 dev(M):dev(M)),
 *     f = q - sigma_y(alpha).
 *
 * The associative flow rule Lp = gamma_dot df/dM is integrated by the
 * exponential update Fp_{n+1} = exp(A) Fp_n with tr(A)=0. The plastic increment
 * A and equivalent plastic increment dgamma are solved simultaneously from the
 * backward-Euler flow and consistency equations. This is a genuine multiplicative
 * finite-strain J2 update, not an additive Green-Lagrange plasticity approximation.
 *
 * All constitutive Jacobians and algorithmic tangents in this implementation are
 * differentiated analytically. The only spectral derivative required by the
 * finite-strain model is the exact Fréchet derivative of exp(A), evaluated from
 * the symmetric eigendecomposition of A through exponential divided differences.
 * No strain or local-unknown finite differences enter the return mapping or the
 * consistent tangent. At hardening-table knots the constitutive law itself is not
 * differentiable; the active tabulated segment supplies the corresponding
 * one-sided tangent.
 */

#include "isotropic_j2_elasticity.h"

#include "../core/logging.h"
#include "strain/axial_strain_green_lagrange.h"
#include "strain/axial_strain_linearized.h"
#include "strain/shell_material_strain_green_lagrange.h"
#include "strain/shell_material_strain_linearized.h"
#include "strain/volume_strain_green_lagrange.h"
#include "strain/volume_strain_linearized.h"
#include "stress/axial_stress_cauchy.h"
#include "stress/axial_stress_pk2.h"
#include "stress/shell_material_stress_cauchy.h"
#include "stress/shell_material_stress_pk2.h"
#include "stress/volume_stress_cauchy.h"
#include "stress/volume_stress_pk2.h"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace fem::material {
namespace {

constexpr Index state_count = 7;
using State      = std::array<Precision, state_count>;
using YieldPoint = IsotropicJ2Elasticity::YieldPoint;
using YieldCurve = std::vector<YieldPoint>;

struct Update3D {
    Mat3  stress = Mat3::Zero();
    State state{};
};

// -----------------------------------------------------------------------------
// Basic tensor/state utilities
// -----------------------------------------------------------------------------

Mat3 sym(const Mat3& A) {
    return Precision(0.5) * (A + A.transpose());
}

Mat3 dev(const Mat3& A) {
    return A - (A.trace() / Precision(3)) * Mat3::Identity();
}

Precision double_contract(const Mat3& A, const Mat3& B) {
    return (A.array() * B.array()).sum();
}

State load_state(const Precision* state) {
    State result{};
    for (Index i = 0; i < state_count; ++i) {
        result[static_cast<std::size_t>(i)] = state[i];
    }
    return result;
}

void store_state(Precision* state, const State& values) {
    for (Index i = 0; i < state_count; ++i) {
        state[i] = values[static_cast<std::size_t>(i)];
    }
}

Mat3 cp_from_state(const State& state) {
    Mat3 Cp;
    Cp << state[0], state[5], state[4],
          state[5], state[1], state[3],
          state[4], state[3], state[2];
    return sym(Cp);
}

void cp_to_state(State& state, const Mat3& Cp) {
    const Mat3 Cps = sym(Cp);
    state[0] = Cps(0, 0);
    state[1] = Cps(1, 1);
    state[2] = Cps(2, 2);
    state[3] = Cps(1, 2);
    state[4] = Cps(0, 2);
    state[5] = Cps(0, 1);
}

Mat3 engineering_strain_tensor(const Vec6& e) {
    Mat3 E;
    E << e(0),                  Precision(0.5) * e(5), Precision(0.5) * e(4),
         Precision(0.5) * e(5), e(1),                  Precision(0.5) * e(3),
         Precision(0.5) * e(4), Precision(0.5) * e(3), e(2);
    return E;
}

Vec6 stress_voigt(const Mat3& S) {
    const Mat3 A = sym(S);
    Vec6 v;
    v << A(0, 0), A(1, 1), A(2, 2), A(1, 2), A(0, 2), A(0, 1);
    return v;
}

template<typename Function>
Mat3 spd_function(const Mat3& A, Function&& function, const char* name) {
    Eigen::SelfAdjointEigenSolver<Mat3> solver(sym(A));
    logging::error(solver.info() == Eigen::Success,
                   "J2: eigendecomposition failed for ", name);

    const Vec3 lambda = solver.eigenvalues();
    const Precision scale = std::max(Precision(1), lambda.cwiseAbs().maxCoeff());
    logging::error(lambda.minCoeff() > Precision(100) * std::numeric_limits<Precision>::epsilon() * scale,
                   "J2: ", name, " is not symmetric positive definite");

    Vec3 mapped;
    for (Index i = 0; i < 3; ++i) {
        mapped(i) = function(lambda(i));
    }

    return sym(solver.eigenvectors() * mapped.asDiagonal() * solver.eigenvectors().transpose());
}

Mat3 spd_sqrt(const Mat3& A) {
    return spd_function(A, [](Precision x) { return std::sqrt(x); }, "metric");
}

Mat3 spd_log(const Mat3& A) {
    return spd_function(A, [](Precision x) { return std::log(x); }, "metric");
}

/**
 * Validates symmetric positive definiteness without evaluating a matrix function.
 * Finite-strain kinematics need this check for C and Ce, but do not need their
 * logarithms; using eigenvalues only avoids the unnecessary spectral mapping.
 */
void require_spd(const Mat3& A, const char* name) {
    Eigen::SelfAdjointEigenSolver<Mat3> solver(sym(A), Eigen::EigenvaluesOnly);
    logging::error(solver.info() == Eigen::Success,
                   "J2: eigendecomposition failed for ", name);

    const Vec3 lambda = solver.eigenvalues();
    const Precision scale = std::max(Precision(1), lambda.cwiseAbs().maxCoeff());
    logging::error(lambda.minCoeff() > Precision(100) * std::numeric_limits<Precision>::epsilon() * scale,
                   "J2: ", name, " is not symmetric positive definite");
}

/**
 * Spectral representation of exp(A) for symmetric A.
 *
 * Keeping the eigenpairs beside exp(A) makes the exact Fréchet derivative cheap:
 * every directional derivative reuses the same eigendecomposition and applies
 * only the exponential divided-difference matrix in the principal basis.
 */
struct SymmetricExponential {
    Mat3 value           = Mat3::Identity();
    Mat3 eigenvectors    = Mat3::Identity();
    Vec3 eigenvalues     = Vec3::Zero();
    Vec3 exp_eigenvalues = Vec3::Ones();
};

SymmetricExponential symmetric_exponential(const Mat3& A) {
    Eigen::SelfAdjointEigenSolver<Mat3> solver(sym(A));
    logging::error(solver.info() == Eigen::Success,
                   "J2: eigendecomposition failed during symmetric exponential");

    SymmetricExponential result;
    result.eigenvalues  = solver.eigenvalues();
    result.eigenvectors = solver.eigenvectors();

    for (Index i = 0; i < 3; ++i) {
        result.exp_eigenvalues(i) = std::exp(result.eigenvalues(i));
    }

    result.value = sym(
        result.eigenvectors
        * result.exp_eigenvalues.asDiagonal()
        * result.eigenvectors.transpose()
    );
    return result;
}

Mat3 sym_exp(const Mat3& A) {
    return symmetric_exponential(A).value;
}

/**
 * Exact Fréchet derivative D exp(A)[H] for symmetric A and H.
 *
 * In the eigenbasis of A the derivative is the Hadamard product
 *
 *     Q^T Dexp(A)[H] Q = L_exp(lambda) .* (Q^T H Q),
 *
 * where the diagonal divided differences are exp(lambda_i) and the off-diagonal
 * entries are (exp(lambda_i)-exp(lambda_j))/(lambda_i-lambda_j). `expm1`
 * evaluates the quotient stably for clustered eigenvalues; the repeated-root
 * limit is exp(lambda).
 */
Mat3 sym_exp_directional(const SymmetricExponential& exponential,
                         const Mat3&                 direction) {
    const Mat3 local_direction =
        exponential.eigenvectors.transpose() * sym(direction) * exponential.eigenvectors;

    Mat3 local_derivative = Mat3::Zero();
    for (Index i = 0; i < 3; ++i) {
        for (Index j = 0; j < 3; ++j) {
            Precision factor = Precision(0);

            if (i == j) {
                factor = exponential.exp_eigenvalues(i);
            } else {
                const Precision lambda_i = exponential.eigenvalues(i);
                const Precision lambda_j = exponential.eigenvalues(j);
                const Precision delta    = lambda_i - lambda_j;
                const Precision scale = std::max({
                    Precision(1), std::abs(lambda_i), std::abs(lambda_j)
                });

                if (std::abs(delta)
                    <= Precision(64) * std::numeric_limits<Precision>::epsilon() * scale) {
                    // Repeated eigenvalue: use the exact divided-difference limit.
                    factor = Precision(0.5)
                           * (exponential.exp_eigenvalues(i)
                              + exponential.exp_eigenvalues(j));
                } else {
                    factor = exponential.exp_eigenvalues(j) * std::expm1(delta) / delta;
                }
            }

            local_derivative(i, j) = factor * local_direction(i, j);
        }
    }

    return sym(
        exponential.eigenvectors
        * local_derivative
        * exponential.eigenvectors.transpose()
    );
}

Mat3 normalize_cp(Mat3 Cp) {
    Cp = sym(Cp);
    const Precision determinant = Cp.determinant();
    logging::error(std::isfinite(determinant) && determinant > Precision(0),
                   "J2: plastic metric lost positive determinant");

    // Associative J2 flow is isochoric. The exponential update preserves det(Cp)=1
    // analytically; remove only accumulated round-off drift.
    Cp /= std::cbrt(determinant);
    return sym(Cp);
}

Precision stress_tolerance(Precision youngs, Precision yield_stress) {
    return Precision(1e-10)
         * std::max({Precision(1), std::abs(youngs), std::abs(yield_stress)});
}

Precision initial_yield_stress(const YieldCurve& yield_curve) {
    logging::error(!yield_curve.empty(),
                   "J2: at least one yield point must be added before evaluation");
    return yield_curve.front().yield_stress;
}

Precision yield_stress_at(const YieldCurve& yield_curve, Precision equivalent_plastic_strain) {
    logging::error(!yield_curve.empty(), "J2: yield curve must not be empty");

    if (yield_curve.size() == 1
        || equivalent_plastic_strain <= yield_curve.front().equivalent_plastic_strain) {
        return yield_curve.front().yield_stress;
    }

    for (std::size_t i = 1; i < yield_curve.size(); ++i) {
        const YieldPoint& left  = yield_curve[i - 1];
        const YieldPoint& right = yield_curve[i];
        if (equivalent_plastic_strain <= right.equivalent_plastic_strain) {
            const Precision span = right.equivalent_plastic_strain
                                 - left.equivalent_plastic_strain;
            const Precision xi = (equivalent_plastic_strain
                                - left.equivalent_plastic_strain) / span;
            return left.yield_stress + xi * (right.yield_stress - left.yield_stress);
        }
    }

    return yield_curve.back().yield_stress;
}

/**
 * Returns the active piecewise-linear hardening slope.
 *
 * At an interior table knot the derivative of the hardening curve is not unique.
 * The interval to the right is selected, which defines the one-sided tangent used
 * by the local Newton system and by the algorithmic material tangent.
 */
Precision hardening_slope_at(const YieldCurve& yield_curve, Precision equivalent_plastic_strain) {
    if (yield_curve.size() < 2
        || equivalent_plastic_strain >= yield_curve.back().equivalent_plastic_strain) {
        return Precision(0);
    }

    for (std::size_t i = 1; i < yield_curve.size(); ++i) {
        const YieldPoint& left  = yield_curve[i - 1];
        const YieldPoint& right = yield_curve[i];
        if (equivalent_plastic_strain < right.equivalent_plastic_strain) {
            return (right.yield_stress - left.yield_stress)
                 / (right.equivalent_plastic_strain - left.equivalent_plastic_strain);
        }
    }

    return Precision(0);
}

Precision solve_small_plastic_increment(Precision q_trial,
                                        Precision alpha_old,
                                        Precision shear,
                                        const YieldCurve& yield_curve) {
    const Precision yield_old = yield_stress_at(yield_curve, alpha_old);
    const Precision f_trial   = q_trial - yield_old;
    if (f_trial <= Precision(0)) {
        return Precision(0);
    }

    Precision lower = Precision(0);
    Precision upper = f_trial / (Precision(3) * shear);
    const Precision tolerance = stress_tolerance(
        Precision(3) * shear,
        initial_yield_stress(yield_curve)
    );

    for (Index iteration = 0; iteration < 80; ++iteration) {
        const Precision gamma = Precision(0.5) * (lower + upper);
        const Precision residual = q_trial
                                 - Precision(3) * shear * gamma
                                 - yield_stress_at(yield_curve, alpha_old + gamma);

        if (std::abs(residual) <= tolerance) {
            return gamma;
        }

        if (residual > Precision(0)) {
            lower = gamma;
        } else {
            upper = gamma;
        }
    }

    return Precision(0.5) * (lower + upper);
}

bool admissible_green_lagrange(const Vec6& e) {
    const Mat3 C = Mat3::Identity() + Precision(2) * engineering_strain_tensor(e);
    Eigen::SelfAdjointEigenSolver<Mat3> solver(sym(C), Eigen::EigenvaluesOnly);
    if (solver.info() != Eigen::Success) {
        return false;
    }
    return solver.eigenvalues().minCoeff() > Precision(1e-12);
}

// -----------------------------------------------------------------------------
// Small-strain three-dimensional constitutive update and exact tangent
// -----------------------------------------------------------------------------

/**
 * Complete small-strain return-map data retained for exact linearization.
 */
struct SmallUpdate {
    Update3D response;
    Mat3     sigma_trial = Mat3::Zero();
    Mat3     s_trial     = Mat3::Zero();
    Mat3     N           = Mat3::Zero();
    Precision q_trial    = Precision(0);
    Precision alpha_old  = Precision(0);
    Precision dgamma     = Precision(0);
    bool      plastic    = false;
};

SmallUpdate integrate_small_strain(const Mat3& strain,
                                   const State& committed,
                                   Precision shear,
                                   Precision bulk,
                                   const YieldCurve& yield_curve) {
    SmallUpdate result;
    result.response.state = committed;

    const Mat3 Cp_old = cp_from_state(committed);
    const Mat3 ep_old = Precision(0.5) * spd_log(Cp_old);
    result.alpha_old  = committed[6];

    const Mat3 ee_trial = sym(strain - ep_old);
    result.sigma_trial = bulk * ee_trial.trace() * Mat3::Identity()
                       + Precision(2) * shear * dev(ee_trial);
    result.s_trial = dev(result.sigma_trial);
    result.q_trial = std::sqrt(
        std::max(
            Precision(0),
            Precision(1.5) * double_contract(result.s_trial, result.s_trial)
        )
    );

    const Precision current_yield = yield_stress_at(yield_curve, result.alpha_old);
    const Precision f_trial = result.q_trial - current_yield;
    result.plastic = f_trial
        > stress_tolerance(Precision(3) * bulk, initial_yield_stress(yield_curve));

    if (!result.plastic) {
        result.response.stress = result.sigma_trial;
        return result;
    }

    logging::error(result.q_trial > Precision(0),
                   "J2: positive yield residual with zero trial equivalent stress");

    result.N = (Precision(1.5) / result.q_trial) * result.s_trial;
    result.dgamma = solve_small_plastic_increment(
        result.q_trial,
        result.alpha_old,
        shear,
        yield_curve
    );

    const Mat3 ep_new = sym(ep_old + result.dgamma * result.N);
    const Mat3 Cp_new = normalize_cp(sym_exp(Precision(2) * ep_new));

    result.response.stress = sym(
        result.sigma_trial - Precision(2) * shear * result.dgamma * result.N
    );
    cp_to_state(result.response.state, Cp_new);
    result.response.state[6] = result.alpha_old + result.dgamma;
    return result;
}

/**
 * Exact small-strain algorithmic tangent of the radial-return update.
 *
 * For a plastic state the scalar consistency equation gives
 *
 *     dgamma = dq_trial / (3 G + H),
 *
 * while the normalized J2 flow direction N = 3/2 s/q is differentiated directly.
 * The six engineering-Voigt columns are assembled from exact tensor directional
 * derivatives; no perturbation of the constitutive update is performed.
 */
Mat6 tangent_small(const Vec6& strain,
                   const SmallUpdate& base,
                   Precision shear,
                   Precision bulk,
                   const YieldCurve& yield_curve) {
    Mat6 C = Mat6::Zero();

    const Precision hardening = base.plastic
        ? hardening_slope_at(yield_curve, base.alpha_old + base.dgamma)
        : Precision(0);

    for (Index column = 0; column < 6; ++column) {
        Vec6 direction = Vec6::Zero();
        direction(column) = Precision(1);

        const Mat3 dstrain = engineering_strain_tensor(direction);
        const Mat3 dsigma_trial = bulk * dstrain.trace() * Mat3::Identity()
                                + Precision(2) * shear * dev(dstrain);

        Mat3 dsigma = dsigma_trial;

        if (base.plastic) {
            const Mat3 ds_trial = dev(dsigma_trial);
            const Precision dq = double_contract(base.N, ds_trial);
            const Precision dgamma = dq / (Precision(3) * shear + hardening);
            const Mat3 dN = (Precision(1.5) / base.q_trial)
                           * (ds_trial - (dq / base.q_trial) * base.s_trial);

            dsigma -= Precision(2) * shear
                    * (dgamma * base.N + base.dgamma * dN);
        }

        C.col(column) = stress_voigt(dsigma);
    }

    return Precision(0.5) * (C + C.transpose());
}

// -----------------------------------------------------------------------------
// Finite-strain three-dimensional constitutive update
// -----------------------------------------------------------------------------

struct FiniteTrial {
    Mat3 A      = Mat3::Zero();
    SymmetricExponential exp_A;
    Mat3 Fp     = Mat3::Identity();
    Mat3 Fp_inv = Mat3::Identity();
    Mat3 Ce     = Mat3::Identity();
    Mat3 Ee     = Mat3::Zero();
    Mat3 Se     = Mat3::Zero();
    Mat3 M      = Mat3::Zero();
    Mat3 m      = Mat3::Zero();
    Mat3 N      = Mat3::Zero();
    Precision q = Precision(0);
};

struct FiniteUpdate {
    Update3D    response;
    Mat3        C         = Mat3::Identity();
    Mat3        Fp_old    = Mat3::Identity();
    Precision   alpha_old = Precision(0);
    Vec6        x         = Vec6::Zero();
    FiniteTrial trial;
    bool        plastic   = false;
};

struct FiniteDirectionalResponse {
    Mat3 stress   = Mat3::Zero();
    Vec6 residual = Vec6::Zero();
};

struct FiniteLocalLinearization {
    Mat6 residual_x = Mat6::Zero();
    Mat6 stress_x   = Mat6::Zero();
};

Mat3 deviatoric_tensor_from_five(const Vec5& v) {
    Mat3 A;
    A << v(0), v(4), v(3),
         v(4), v(1), v(2),
         v(3), v(2), -v(0) - v(1);
    return sym(A);
}

Vec5 five_from_deviatoric_tensor(const Mat3& A) {
    const Mat3 D = dev(sym(A));
    Vec5 v;
    v << D(0, 0), D(1, 1), D(1, 2), D(0, 2), D(0, 1);
    return v;
}

Mat3 finite_a_direction(Index component) {
    Vec5 direction = Vec5::Zero();
    direction(component) = Precision(1);
    return deviatoric_tensor_from_five(direction);
}

FiniteTrial finite_trial_from_increment(const Mat3& C,
                                        const Mat3& Fp_old,
                                        const Vec6& x,
                                        Precision shear,
                                        Precision bulk) {
    FiniteTrial trial;

    trial.A     = deviatoric_tensor_from_five(x.template head<5>());
    trial.exp_A = symmetric_exponential(trial.A);

    trial.Fp     = trial.exp_A.value * Fp_old;
    trial.Fp_inv = trial.Fp.inverse();
    trial.Ce     = sym(trial.Fp_inv.transpose() * C * trial.Fp_inv);
    require_spd(trial.Ce, "elastic metric");

    trial.Ee = Precision(0.5) * (trial.Ce - Mat3::Identity());
    const Precision lame = bulk - Precision(2.0 / 3.0) * shear;
    trial.Se = lame * trial.Ee.trace() * Mat3::Identity()
             + Precision(2) * shear * trial.Ee;

    trial.M = sym(trial.Ce * trial.Se);
    trial.m = dev(trial.M);
    trial.q = std::sqrt(
        std::max(Precision(0), Precision(1.5) * double_contract(trial.m, trial.m))
    );

    if (trial.q > Precision(0)) {
        trial.N = (Precision(1.5) / trial.q) * trial.m;
    }

    return trial;
}

Mat3 finite_pk2_stress(const FiniteTrial& trial) {
    return sym(trial.Fp_inv * trial.Se * trial.Fp_inv.transpose());
}

Vec6 finite_local_residual(const FiniteTrial& trial,
                           Precision alpha_old,
                           const Vec6& x,
                           Precision shear,
                           const YieldCurve& yield_curve) {
    const Precision dgamma = x(5);

    logging::error(dgamma >= Precision(-1e-12),
                   "J2: negative equivalent plastic increment in local solve");
    logging::error(trial.q > Precision(0),
                   "J2: plastic local solve reached zero Mandel equivalent stress");

    Vec6 residual = Vec6::Zero();
    residual.template head<5>() = five_from_deviatoric_tensor(
        trial.A - dgamma * trial.N
    );

    const Precision alpha = alpha_old + dgamma;
    const Precision hardening = hardening_slope_at(yield_curve, alpha);
    const Precision denom = Precision(3) * shear + hardening;
    residual(5) = (trial.q - yield_stress_at(yield_curve, alpha)) / denom;
    return residual;
}

Vec6 finite_local_residual(const Mat3& C,
                           const Mat3& Fp_old,
                           Precision alpha_old,
                           const Vec6& x,
                           Precision shear,
                           Precision bulk,
                           const YieldCurve& yield_curve) {
    const FiniteTrial trial = finite_trial_from_increment(C, Fp_old, x, shear, bulk);
    return finite_local_residual(trial, alpha_old, x, shear, yield_curve);
}

/**
 * Exact directional derivative of finite-strain stress and, optionally, of the
 * six local return-mapping equations.
 *
 * The supplied directions represent independent variations dC, dA and dgamma.
 * The exponential-map contribution is differentiated by the exact Fréchet
 * derivative Dexp(A)[dA]. All remaining quantities follow by the ordinary chain
 * rule:
 *
 *     dFp^{-1} = -Fp^{-1} dFp Fp^{-1},
 *     dCe      = sym(dFp^{-T} C Fp^{-1}
 *                    + Fp^{-T} dC Fp^{-1}
 *                    + Fp^{-T} C dFp^{-1}),
 *     dq       = N : dM,
 *     dN       = 3/(2q) [dev(dM) - (dq/q) dev(M)].
 *
 * This function is the common analytic building block for both the local Newton
 * Jacobian and the global consistent material tangent.
 */
FiniteDirectionalResponse finite_directional_response(
    const FiniteTrial& trial,
    const Mat3&        C,
    const Mat3&        Fp_old,
    Precision          alpha_old,
    const Vec6&        x,
    const Mat3&        dC,
    const Mat3&        dA,
    Precision          dgamma,
    Precision          shear,
    Precision          bulk,
    const YieldCurve&  yield_curve,
    bool               include_residual
) {
    const Mat3 dplastic_increment = sym_exp_directional(trial.exp_A, dA);
    const Mat3 dFp = dplastic_increment * Fp_old;
    const Mat3 dFp_inv = -trial.Fp_inv * dFp * trial.Fp_inv;

    // Differentiate Ce from its defining expression. Keeping the original total
    // metric C explicit avoids algebraic cancellation assumptions and keeps the
    // derivative valid for simultaneous dC and dA variations.
    const Mat3 dCe = sym(
        dFp_inv.transpose() * C * trial.Fp_inv
        + trial.Fp_inv.transpose() * dC * trial.Fp_inv
        + trial.Fp_inv.transpose() * C * dFp_inv
    );

    const Mat3 dEe = Precision(0.5) * dCe;
    const Precision lame = bulk - Precision(2.0 / 3.0) * shear;
    const Mat3 dSe = lame * dEe.trace() * Mat3::Identity()
                   + Precision(2) * shear * dEe;

    const Mat3 dM = sym(dCe * trial.Se + trial.Ce * dSe);
    const Mat3 dm = dev(dM);

    Precision dq = Precision(0);
    Mat3 dN = Mat3::Zero();
    if (trial.q > Precision(0)) {
        dq = double_contract(trial.N, dm);
        dN = (Precision(1.5) / trial.q)
           * (dm - (dq / trial.q) * trial.m);
    }

    FiniteDirectionalResponse result;
    result.stress = sym(
        dFp_inv * trial.Se * trial.Fp_inv.transpose()
        + trial.Fp_inv * dSe * trial.Fp_inv.transpose()
        + trial.Fp_inv * trial.Se * dFp_inv.transpose()
    );

    if (!include_residual) {
        return result;
    }

    logging::error(trial.q > Precision(0),
                   "J2: cannot differentiate plastic residual at zero Mandel stress");

    result.residual.template head<5>() = five_from_deviatoric_tensor(
        dA - dgamma * trial.N - x(5) * dN
    );

    const Precision alpha = alpha_old + x(5);
    const Precision hardening = hardening_slope_at(yield_curve, alpha);
    const Precision denom = Precision(3) * shear + hardening;

    // H is constant inside each tabulated segment. At a knot the active right
    // segment defines the one-sided derivative by hardening_slope_at().
    result.residual(5) = (dq - hardening * dgamma) / denom;
    return result;
}

/**
 * Builds the exact local Jacobian dR/dx and, at the same time, dS/dx.
 *
 * The five tensor unknowns use the exact linear basis of the symmetric
 * trace-free A parameterization. The sixth column is the direct derivative with
 * respect to dgamma. No perturbed local solves or finite-difference probes are
 * required.
 */
FiniteLocalLinearization finite_local_linearization_x(
    const FiniteTrial& trial,
    const Mat3&        C,
    const Mat3&        Fp_old,
    Precision          alpha_old,
    const Vec6&        x,
    Precision          shear,
    Precision          bulk,
    const YieldCurve&  yield_curve
) {
    FiniteLocalLinearization result;

    for (Index column = 0; column < 6; ++column) {
        const Mat3 dA = column < 5
            ? finite_a_direction(column)
            : Mat3::Zero();
        const Precision dgamma = column == 5 ? Precision(1) : Precision(0);

        const FiniteDirectionalResponse directional = finite_directional_response(
            trial,
            C,
            Fp_old,
            alpha_old,
            x,
            Mat3::Zero(),
            dA,
            dgamma,
            shear,
            bulk,
            yield_curve,
            true
        );

        result.residual_x.col(column) = directional.residual;
        result.stress_x.col(column)   = stress_voigt(directional.stress);
    }

    return result;
}

FiniteUpdate integrate_finite_strain(const Mat3& green_lagrange,
                                     const State& committed,
                                     Precision shear,
                                     Precision bulk,
                                     const YieldCurve& yield_curve) {
    FiniteUpdate result;
    result.response.state = committed;

    result.C = sym(Mat3::Identity() + Precision(2) * green_lagrange);
    require_spd(result.C, "total metric");

    const Mat3 Cp_old = cp_from_state(committed);
    result.Fp_old    = spd_sqrt(Cp_old);
    result.alpha_old = committed[6];
    result.x.setZero();

    result.trial = finite_trial_from_increment(
        result.C, result.Fp_old, result.x, shear, bulk
    );

    const Precision current_yield = yield_stress_at(yield_curve, result.alpha_old);
    const Precision f_trial = result.trial.q - current_yield;
    result.plastic = f_trial > stress_tolerance(
        Precision(3) * bulk,
        initial_yield_stress(yield_curve)
    );

    if (result.plastic) {
        logging::error(result.trial.q > Precision(0),
                       "J2: positive yield residual with zero trial Mandel stress");

        // The elastic trial is the natural starting state. The local Newton
        // Jacobian below is the exact derivative of the six coupled backward-
        // Euler equations with respect to A and dgamma.
        result.x.setZero();

        bool local_converged = false;
        for (Index iteration = 0; iteration < 30; ++iteration) {
            const FiniteTrial iterate_trial = finite_trial_from_increment(
                result.C, result.Fp_old, result.x, shear, bulk
            );
            const Vec6 residual = finite_local_residual(
                iterate_trial, result.alpha_old, result.x, shear, yield_curve
            );

            if (residual.cwiseAbs().maxCoeff() <= Precision(1e-10)) {
                result.trial = iterate_trial;
                local_converged = true;
                break;
            }

            const FiniteLocalLinearization linearization = finite_local_linearization_x(
                iterate_trial,
                result.C,
                result.Fp_old,
                result.alpha_old,
                result.x,
                shear,
                bulk,
                yield_curve
            );

            Eigen::FullPivLU<Mat6> lu(linearization.residual_x);
            logging::error(lu.isInvertible(),
                           "J2: singular local finite-plasticity Jacobian");
            const Vec6 delta = lu.solve(-residual);

            Precision scale = Precision(1);
            bool accepted = false;
            const Precision residual_norm = residual.norm();
            while (scale >= Precision(1e-8)) {
                const Vec6 candidate = result.x + scale * delta;
                if (candidate(5) >= Precision(0)) {
                    const Vec6 candidate_residual = finite_local_residual(
                        result.C,
                        result.Fp_old,
                        result.alpha_old,
                        candidate,
                        shear,
                        bulk,
                        yield_curve
                    );
                    if (candidate_residual.norm() < residual_norm) {
                        result.x = candidate;
                        accepted = true;
                        break;
                    }
                }
                scale *= Precision(0.5);
            }

            logging::error(accepted,
                           "J2: finite-plasticity local Newton line search failed");
        }

        logging::error(local_converged,
                       "J2: finite-plasticity local return mapping did not converge");
    }

    result.response.stress = finite_pk2_stress(result.trial);

    const Mat3 Cp_new = normalize_cp(result.trial.Fp.transpose() * result.trial.Fp);
    cp_to_state(result.response.state, Cp_new);
    result.response.state[6] = result.alpha_old + result.x(5);
    return result;
}

// -----------------------------------------------------------------------------
// Exact finite-strain algorithmic tangent
// -----------------------------------------------------------------------------

/**
 * Builds the exact consistent tangent of the converged finite return map.
 *
 * For plastic states the local equations R(x,E)=0 give
 *
 *     dx/dE = -(dR/dx)^-1 dR/dE,
 *
 * and therefore
 *
 *     dS/dE = dS/dE|x + dS/dx dx/dE.
 *
 * Every partial derivative on the right-hand side is obtained by the analytic
 * tensor directional derivative above. In particular Dexp(A) is evaluated by its
 * exact spectral Fréchet derivative. Elastic states require only dS/dE at frozen
 * plastic history.
 */
Mat6 tangent_finite(const Vec6& strain,
                    const FiniteUpdate& base,
                    Precision shear,
                    Precision bulk,
                    const YieldCurve& yield_curve) {
    Mat6 stress_E   = Mat6::Zero();
    Mat6 residual_E = Mat6::Zero();

    // C = I + 2E. Each engineering-Voigt strain basis therefore produces the
    // exact tensor direction dC = 2 dE.
    for (Index column = 0; column < 6; ++column) {
        Vec6 direction = Vec6::Zero();
        direction(column) = Precision(1);
        const Mat3 dC = Precision(2) * engineering_strain_tensor(direction);

        const FiniteDirectionalResponse directional = finite_directional_response(
            base.trial,
            base.C,
            base.Fp_old,
            base.alpha_old,
            base.x,
            dC,
            Mat3::Zero(),
            Precision(0),
            shear,
            bulk,
            yield_curve,
            base.plastic
        );

        stress_E.col(column) = stress_voigt(directional.stress);
        if (base.plastic) {
            residual_E.col(column) = directional.residual;
        }
    }

    if (!base.plastic) {
        return Precision(0.5) * (stress_E + stress_E.transpose());
    }

    const FiniteLocalLinearization local = finite_local_linearization_x(
        base.trial,
        base.C,
        base.Fp_old,
        base.alpha_old,
        base.x,
        shear,
        bulk,
        yield_curve
    );

    Eigen::FullPivLU<Mat6> lu(local.residual_x);
    logging::error(lu.isInvertible(),
                   "J2: singular local Jacobian during finite tangent linearization");

    const Mat6 dx_dE = lu.solve(-residual_E);
    const Mat6 C = stress_E + local.stress_x * dx_dE;

    // Associative J2 plasticity with isotropic elasticity has a symmetric
    // algorithmic tangent. Remove only round-off asymmetry accumulated in the
    // spectral and dense local algebra.
    return Precision(0.5) * (C + C.transpose());
}

// -----------------------------------------------------------------------------
// Plane-stress and uniaxial-stress reductions
// -----------------------------------------------------------------------------

Vec6 shell_to_volume_strain(const Vec5& shell, Precision e33) {
    Vec6 e;
    e << shell(0), shell(1), e33, shell(4), shell(3), shell(2);
    return e;
}

Vec5 volume_to_shell_stress(const Vec6& volume) {
    Vec5 shell;
    shell << volume(0), volume(1), volume(5), volume(4), volume(3);
    return shell;
}

Mat5 condense_shell_tangent(const Mat6& C) {
    constexpr std::array<Index, 5> a{0, 1, 5, 4, 3};
    constexpr Index z = 2;

    Mat5 Caa;
    Eigen::Matrix<Precision, 5, 1> Caz;
    Eigen::Matrix<Precision, 1, 5> Cza;

    for (Index i = 0; i < 5; ++i) {
        Caz(i) = C(a[static_cast<std::size_t>(i)], z);
        Cza(i) = C(z, a[static_cast<std::size_t>(i)]);
        for (Index j = 0; j < 5; ++j) {
            Caa(i, j) = C(a[static_cast<std::size_t>(i)], a[static_cast<std::size_t>(j)]);
        }
    }

    logging::error(std::abs(C(z, z)) > Precision(100) * std::numeric_limits<Precision>::epsilon(),
                   "J2: singular thickness-normal tangent during shell condensation");

    return Caa - (Caz * Cza) / C(z, z);
}

Precision condense_axial_tangent(const Mat6& C) {
    Eigen::Matrix<Precision, 1, 2> Cab;
    Eigen::Matrix<Precision, 2, 1> Cba;
    Eigen::Matrix<Precision, 2, 2> Cbb;

    Cab << C(0, 1), C(0, 2);
    Cba << C(1, 0), C(2, 0);
    Cbb << C(1, 1), C(1, 2),
           C(2, 1), C(2, 2);

    Eigen::FullPivLU<Eigen::Matrix<Precision, 2, 2>> lu(Cbb);
    logging::error(lu.isInvertible(),
                   "J2: singular transverse block during axial tangent condensation");

    return C(0, 0) - (Cab * lu.solve(Cba))(0, 0);
}

/**
 * Evaluates one full three-dimensional constitutive candidate together with its
 * exact algorithmic tangent from the same committed state.
 */
template<bool Finite>
Update3D evaluate_volume_candidate(const Vec6& strain,
                                   const State& committed,
                                   Precision shear,
                                   Precision bulk,
                                   const YieldCurve& yield_curve,
                                   Mat6& tangent) {
    if constexpr (Finite) {
        const FiniteUpdate update = integrate_finite_strain(
            engineering_strain_tensor(strain), committed,
            shear, bulk, yield_curve
        );
        tangent = tangent_finite(strain, update, shear, bulk, yield_curve);
        return update.response;
    } else {
        const SmallUpdate update = integrate_small_strain(
            engineering_strain_tensor(strain), committed,
            shear, bulk, yield_curve
        );
        tangent = tangent_small(strain, update, shear, bulk, yield_curve);
        return update.response;
    }
}

/**
 * Solves S33 = 0 with the exact three-dimensional consistent tangent.
 *
 * The Newton derivative is simply C3333 of the same constitutive candidate; no
 * perturbation of e33 and no auxiliary return maps are required.
 */
template<bool Finite>
Update3D solve_shell_plane_stress(const Vec5& shell_strain,
                                  const State& committed,
                                  Precision youngs,
                                  Precision poisson,
                                  Precision shear,
                                  Precision bulk,
                                  const YieldCurve& yield_curve,
                                  Vec6& converged_strain,
                                  Mat6& converged_tangent) {
    Precision e33 = -poisson / (Precision(1) - poisson)
                  * (shell_strain(0) + shell_strain(1));

    if constexpr (Finite) {
        e33 = std::max(e33, Precision(-0.49));
    }

    const Precision tolerance = stress_tolerance(youngs, initial_yield_stress(yield_curve));
    Update3D current;
    bool converged = false;

    for (Index iteration = 0; iteration < 30; ++iteration) {
        const Vec6 e = shell_to_volume_strain(shell_strain, e33);
        if constexpr (Finite) {
            logging::error(admissible_green_lagrange(e),
                           "J2: shell plane-stress iterate left admissible Green-Lagrange domain");
        }

        Mat6 current_tangent;
        current = evaluate_volume_candidate<Finite>(
            e, committed, shear, bulk, yield_curve, current_tangent
        );

        const Precision residual = current.stress(2, 2);
        if (std::abs(residual) <= tolerance) {
            converged = true;
            converged_strain  = e;
            converged_tangent = current_tangent;
            break;
        }

        const Precision derivative = current_tangent(2, 2);
        logging::error(std::isfinite(derivative)
                       && std::abs(derivative) > Precision(100) * std::numeric_limits<Precision>::epsilon(),
                       "J2: singular shell plane-stress Newton derivative");

        const Precision delta = -residual / derivative;

        if constexpr (Finite) {
            Precision scale = Precision(1);
            bool accepted = false;
            while (scale >= Precision(1e-8)) {
                const Vec6 candidate = shell_to_volume_strain(shell_strain, e33 + scale * delta);
                if (admissible_green_lagrange(candidate)) {
                    e33 += scale * delta;
                    accepted = true;
                    break;
                }
                scale *= Precision(0.5);
            }
            logging::error(accepted,
                           "J2: shell plane-stress Newton could not find admissible step");
        } else {
            e33 += delta;
        }
    }

    logging::error(converged,
                   "J2: shell plane-stress constitutive solve did not converge");
    return current;
}

/**
 * Solves S22 = S33 = 0 with the exact transverse block of the consistent
 * three-dimensional tangent.
 */
template<bool Finite>
Update3D solve_axial_stress(const Precision axial_strain,
                            const State& committed,
                            Precision youngs,
                            Precision poisson,
                            Precision shear,
                            Precision bulk,
                            const YieldCurve& yield_curve,
                            Vec6& converged_strain,
                            Mat6& converged_tangent) {
    Vec6 e = Vec6::Zero();
    e(0) = axial_strain;
    e(1) = -poisson * axial_strain;
    e(2) = -poisson * axial_strain;

    if constexpr (Finite) {
        e(1) = std::max(e(1), Precision(-0.49));
        e(2) = std::max(e(2), Precision(-0.49));
    }

    const Precision tolerance = stress_tolerance(youngs, initial_yield_stress(yield_curve));
    Update3D current;
    bool converged = false;

    for (Index iteration = 0; iteration < 30; ++iteration) {
        if constexpr (Finite) {
            logging::error(admissible_green_lagrange(e),
                           "J2: axial constitutive iterate left admissible Green-Lagrange domain");
        }

        Mat6 current_tangent;
        current = evaluate_volume_candidate<Finite>(
            e, committed, shear, bulk, yield_curve, current_tangent
        );

        Eigen::Matrix<Precision, 2, 1> residual;
        residual << current.stress(1, 1), current.stress(2, 2);

        if (residual.cwiseAbs().maxCoeff() <= tolerance) {
            converged = true;
            converged_strain  = e;
            converged_tangent = current_tangent;
            break;
        }

        Eigen::Matrix<Precision, 2, 2> J;
        J << current_tangent(1, 1), current_tangent(1, 2),
             current_tangent(2, 1), current_tangent(2, 2);

        Eigen::FullPivLU<Eigen::Matrix<Precision, 2, 2>> lu(J);
        logging::error(lu.isInvertible(),
                       "J2: singular transverse Jacobian in axial stress reduction");
        const Eigen::Matrix<Precision, 2, 1> delta = lu.solve(-residual);

        if constexpr (Finite) {
            Precision scale = Precision(1);
            bool accepted = false;
            while (scale >= Precision(1e-8)) {
                Vec6 candidate = e;
                candidate(1) += scale * delta(0);
                candidate(2) += scale * delta(1);
                if (admissible_green_lagrange(candidate)) {
                    e = candidate;
                    accepted = true;
                    break;
                }
                scale *= Precision(0.5);
            }
            logging::error(accepted,
                           "J2: axial constitutive Newton could not find admissible step");
        } else {
            e(1) += delta(0);
            e(2) += delta(1);
        }
    }

    logging::error(converged,
                   "J2: axial transverse-stress constitutive solve did not converge");
    return current;
}

} // namespace

// -----------------------------------------------------------------------------
// Public material implementation
// -----------------------------------------------------------------------------

IsotropicJ2Elasticity::IsotropicJ2Elasticity(Precision youngs_in, Precision poisson_in)
    : youngs  (youngs_in),
      poisson (poisson_in) {
    logging::error(youngs > Precision(0),
                   "J2: Young's modulus must be positive");
    logging::error(poisson > Precision(-1) && poisson < Precision(0.5),
                   "J2: Poisson ratio must be in (-1, 0.5)");
}

Precision IsotropicJ2Elasticity::shear_modulus() const {
    return youngs / (Precision(2) * (Precision(1) + poisson));
}

Precision IsotropicJ2Elasticity::bulk_modulus() const {
    return youngs / (Precision(3) * (Precision(1) - Precision(2) * poisson));
}

void IsotropicJ2Elasticity::add_yield_point(Precision yield_stress,
                                            Precision equivalent_plastic_strain) {
    logging::error(std::isfinite(yield_stress) && yield_stress > Precision(0),
                   "J2: yield stress must be finite and positive");
    logging::error(std::isfinite(equivalent_plastic_strain)
                   && equivalent_plastic_strain >= Precision(0),
                   "J2: equivalent plastic strain must be finite and non-negative");

    const Precision strain_tolerance = Precision(100)
        * std::numeric_limits<Precision>::epsilon();

    if (yield_points_.empty()) {
        logging::error(std::abs(equivalent_plastic_strain) <= strain_tolerance,
                       "J2: first yield point must have zero equivalent plastic strain");
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

const std::vector<IsotropicJ2Elasticity::YieldPoint>&
IsotropicJ2Elasticity::get_yield_points() const {
    return yield_points_;
}

bool IsotropicJ2Elasticity::supports_axial_linearized() const { return true; }
bool IsotropicJ2Elasticity::supports_axial_green_lagrange() const { return true; }
bool IsotropicJ2Elasticity::supports_volume_linearized() const { return true; }
bool IsotropicJ2Elasticity::supports_volume_green_lagrange() const { return true; }
bool IsotropicJ2Elasticity::supports_shell_integration_linearized() const { return true; }
bool IsotropicJ2Elasticity::supports_shell_integration_green_lagrange() const { return true; }

Index IsotropicJ2Elasticity::state_size() const {
    return state_count;
}

void IsotropicJ2Elasticity::initialize_state(Precision* state) const {
    state[0] = Precision(1);
    state[1] = Precision(1);
    state[2] = Precision(1);
    state[3] = Precision(0);
    state[4] = Precision(0);
    state[5] = Precision(0);
    state[6] = Precision(0);
}

void IsotropicJ2Elasticity::evaluate(const VolumeStrainLinearized& strain,
                                     Precision*                    state,
                                     VolumeStressCauchy&           stress,
                                     Mat6&                         tangent) const {
    const State committed = load_state(state);
    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();
    const SmallUpdate update = integrate_small_strain(
        strain.tensor(), committed, shear, bulk, yield_points_
    );

    stress  = VolumeStressCauchy(update.response.stress);
    tangent = tangent_small(strain.voigt(), update, shear, bulk, yield_points_);
    store_state(state, update.response.state);
}

void IsotropicJ2Elasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                     Precision*                       state,
                                     VolumeStressPK2&                 stress,
                                     Mat6*                            tangent) const {
    const State committed = load_state(state);
    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();
    const FiniteUpdate update = integrate_finite_strain(
        strain.tensor(), committed, shear, bulk, yield_points_
    );

    stress = VolumeStressPK2(update.response.stress);
    if (tangent != nullptr) {
        *tangent = tangent_finite(
            strain.voigt(), update, shear, bulk, yield_points_
        );
    }
    store_state(state, update.response.state);
}

void IsotropicJ2Elasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                     Precision*                           state,
                                     ShellMaterialStressCauchy&           stress,
                                     Mat5&                                tangent) const {
    const State committed = load_state(state);
    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();

    Vec6 volume_strain;
    Mat6 C3;
    const Update3D update = solve_shell_plane_stress<false>(
        strain.values(),
        committed,
        youngs,
        poisson,
        shear,
        bulk,
        yield_points_,
        volume_strain,
        C3
    );

    stress.values() = volume_to_shell_stress(stress_voigt(update.stress));
    tangent = condense_shell_tangent(C3);
    store_state(state, update.state);
}

void IsotropicJ2Elasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                                     Precision*                              state,
                                     ShellMaterialStressPK2&                 stress,
                                     Mat5&                                   tangent) const {
    const State committed = load_state(state);
    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();

    Vec6 volume_strain;
    Mat6 C3;
    const Update3D update = solve_shell_plane_stress<true>(
        strain.values(),
        committed,
        youngs,
        poisson,
        shear,
        bulk,
        yield_points_,
        volume_strain,
        C3
    );

    stress.values() = volume_to_shell_stress(stress_voigt(update.stress));
    tangent = condense_shell_tangent(C3);
    store_state(state, update.state);
}

void IsotropicJ2Elasticity::evaluate(const AxialStrainLinearized& strain,
                                     Precision*                   state,
                                     AxialStressCauchy&           stress,
                                     Precision&                   tangent) const {
    const State committed = load_state(state);
    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();

    Vec6 volume_strain;
    Mat6 C3;
    const Update3D update = solve_axial_stress<false>(
        strain.value(),
        committed,
        youngs,
        poisson,
        shear,
        bulk,
        yield_points_,
        volume_strain,
        C3
    );

    stress.value() = update.stress(0, 0);
    tangent = condense_axial_tangent(C3);
    store_state(state, update.state);
}

void IsotropicJ2Elasticity::evaluate(const AxialStrainGreenLagrange& strain,
                                     Precision*                      state,
                                     AxialStressPK2&                 stress,
                                     Precision&                      tangent) const {
    const State committed = load_state(state);
    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();

    Vec6 volume_strain;
    Mat6 C3;
    const Update3D update = solve_axial_stress<true>(
        strain.value(),
        committed,
        youngs,
        poisson,
        shear,
        bulk,
        yield_points_,
        volume_strain,
        C3
    );

    stress.value() = update.stress(0, 0);
    tangent = condense_axial_tangent(C3);
    store_state(state, update.state);
}

} // namespace fem::material
