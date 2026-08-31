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
 * The returned material tangents are numerical derivatives of the complete
 * converged discrete constitutive update from the same committed state. The
 * stress/history update itself is the fully implicit finite-strain model; a later
 * closed-form linearization can replace the numerical derivative without changing
 * the constitutive equations.
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
using State = std::array<Precision, state_count>;
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

Mat3 sym_exp(const Mat3& A) {
    Eigen::SelfAdjointEigenSolver<Mat3> solver(sym(A));
    logging::error(solver.info() == Eigen::Success,
                   "J2: eigendecomposition failed during symmetric exponential");

    Vec3 mapped;
    for (Index i = 0; i < 3; ++i) {
        mapped(i) = std::exp(solver.eigenvalues()(i));
    }

    return sym(solver.eigenvectors() * mapped.asDiagonal() * solver.eigenvectors().transpose());
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

Precision derivative_step(Precision value) {
    return Precision(1e-7) * std::max(Precision(1), std::abs(value));
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
// Three-dimensional constitutive kernels
// -----------------------------------------------------------------------------

Update3D update_small_strain(const Mat3& strain,
                             const State& committed,
                             Precision shear,
                             Precision bulk,
                             const YieldCurve& yield_curve) {
    Update3D result;
    result.state = committed;

    const Mat3 Cp_old = cp_from_state(committed);
    const Mat3 ep_old = Precision(0.5) * spd_log(Cp_old);
    const Precision alpha_old = committed[6];

    const Mat3 ee_trial = sym(strain - ep_old);
    const Mat3 sigma_trial = bulk * ee_trial.trace() * Mat3::Identity()
                           + Precision(2) * shear * dev(ee_trial);
    const Mat3 s_trial = dev(sigma_trial);
    const Precision q_trial = std::sqrt(
        std::max(Precision(0), Precision(1.5) * double_contract(s_trial, s_trial))
    );
    const Precision current_yield = yield_stress_at(yield_curve, alpha_old);
    const Precision f_trial = q_trial - current_yield;

    if (f_trial <= stress_tolerance(Precision(3) * bulk, initial_yield_stress(yield_curve))) {
        result.stress = sigma_trial;
        return result;
    }

    logging::error(q_trial > Precision(0),
                   "J2: positive yield residual with zero trial equivalent stress");

    const Mat3 N = (Precision(1.5) / q_trial) * s_trial;
    const Precision dgamma = solve_small_plastic_increment(
        q_trial, alpha_old, shear, yield_curve
    );

    const Mat3 ep_new = sym(ep_old + dgamma * N);
    const Mat3 Cp_new = normalize_cp(sym_exp(Precision(2) * ep_new));

    result.stress = sym(sigma_trial - Precision(2) * shear * dgamma * N);
    cp_to_state(result.state, Cp_new);
    result.state[6] = alpha_old + dgamma;
    return result;
}

struct FiniteTrial {
    Mat3 Fp = Mat3::Identity();
    Mat3 Ce = Mat3::Identity();
    Mat3 Ee = Mat3::Zero();
    Mat3 Se = Mat3::Zero();
    Mat3 M  = Mat3::Zero();
    Precision q = Precision(0);
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

FiniteTrial finite_trial_from_increment(const Mat3& C,
                                        const Mat3& Fp_old,
                                        const Vec6& x,
                                        Precision shear,
                                        Precision bulk) {
    FiniteTrial trial;

    const Vec5 a = x.template head<5>();
    const Mat3 A = deviatoric_tensor_from_five(a);
    const Mat3 plastic_increment = sym_exp(A);

    trial.Fp = plastic_increment * Fp_old;
    const Mat3 Fp_inv = trial.Fp.inverse();
    trial.Ce = sym(Fp_inv.transpose() * C * Fp_inv);
    (void) spd_log(trial.Ce);

    trial.Ee = Precision(0.5) * (trial.Ce - Mat3::Identity());
    const Precision lame = bulk - Precision(2.0 / 3.0) * shear;
    trial.Se = lame * trial.Ee.trace() * Mat3::Identity()
             + Precision(2) * shear * trial.Ee;

    trial.M = sym(trial.Ce * trial.Se);
    const Mat3 m = dev(trial.M);
    trial.q = std::sqrt(
        std::max(Precision(0), Precision(1.5) * double_contract(m, m))
    );
    return trial;
}

Vec6 finite_local_residual(const Mat3& C,
                           const Mat3& Fp_old,
                           Precision alpha_old,
                           const Vec6& x,
                           Precision shear,
                           Precision bulk,
                           const YieldCurve& yield_curve) {
    const FiniteTrial trial = finite_trial_from_increment(C, Fp_old, x, shear, bulk);
    const Precision dgamma = x(5);

    Vec6 residual = Vec6::Zero();
    const Mat3 A = deviatoric_tensor_from_five(x.template head<5>());

    logging::error(dgamma >= Precision(-1e-12),
                   "J2: negative equivalent plastic increment in local solve");
    logging::error(trial.q > Precision(0),
                   "J2: plastic local solve reached zero Mandel equivalent stress");

    const Mat3 N = (Precision(1.5) / trial.q) * dev(trial.M);
    residual.template head<5>() = five_from_deviatoric_tensor(A - dgamma * N);

    const Precision alpha = alpha_old + dgamma;
    const Precision local_hardening = hardening_slope_at(yield_curve, alpha);
    const Precision denom = Precision(3) * shear + local_hardening;
    residual(5) = (trial.q - yield_stress_at(yield_curve, alpha)) / denom;
    return residual;
}

Update3D update_finite_strain(const Mat3& green_lagrange,
                              const State& committed,
                              Precision shear,
                              Precision bulk,
                              const YieldCurve& yield_curve) {
    Update3D result;
    result.state = committed;

    const Mat3 C = sym(Mat3::Identity() + Precision(2) * green_lagrange);
    (void) spd_log(C);

    const Mat3 Cp_old = cp_from_state(committed);
    const Mat3 Fp_old = spd_sqrt(Cp_old);
    const Precision alpha_old = committed[6];

    Vec6 x = Vec6::Zero();
    const FiniteTrial elastic_trial = finite_trial_from_increment(
        C, Fp_old, x, shear, bulk
    );
    const Precision current_yield = yield_stress_at(yield_curve, alpha_old);
    const Precision f_trial = elastic_trial.q - current_yield;

    FiniteTrial converged = elastic_trial;

    if (f_trial > stress_tolerance(Precision(3) * bulk, initial_yield_stress(yield_curve))) {
        logging::error(elastic_trial.q > Precision(0),
                       "J2: positive yield residual with zero trial Mandel stress");

        const Mat3 N_trial = (Precision(1.5) / elastic_trial.q) * dev(elastic_trial.M);
        const Precision gamma0 = f_trial
            / (Precision(3) * shear + hardening_slope_at(yield_curve, alpha_old));
        x.template head<5>() = five_from_deviatoric_tensor(gamma0 * N_trial);
        x(5) = gamma0;

        bool local_converged = false;
        for (Index iteration = 0; iteration < 30; ++iteration) {
            const Vec6 residual = finite_local_residual(
                C, Fp_old, alpha_old, x,
                shear, bulk, yield_curve
            );

            if (residual.cwiseAbs().maxCoeff() <= Precision(1e-10)) {
                local_converged = true;
                break;
            }

            Mat6 J = Mat6::Zero();
            for (Index column = 0; column < 6; ++column) {
                const Precision h = derivative_step(x(column));
                Vec6 xp = x;
                Vec6 xm = x;
                xp(column) += h;
                xm(column) -= h;

                if (column == 5 && xm(5) < Precision(0)) {
                    const Vec6 rp = finite_local_residual(
                        C, Fp_old, alpha_old, xp,
                        shear, bulk, yield_curve
                    );
                    J.col(column) = (rp - residual) / h;
                } else {
                    const Vec6 rp = finite_local_residual(
                        C, Fp_old, alpha_old, xp,
                        shear, bulk, yield_curve
                    );
                    const Vec6 rm = finite_local_residual(
                        C, Fp_old, alpha_old, xm,
                        shear, bulk, yield_curve
                    );
                    J.col(column) = (rp - rm) / (Precision(2) * h);
                }
            }

            Eigen::FullPivLU<Mat6> lu(J);
            logging::error(lu.isInvertible(),
                           "J2: singular local finite-plasticity Jacobian");
            const Vec6 delta = lu.solve(-residual);

            Precision scale = Precision(1);
            bool accepted = false;
            const Precision residual_norm = residual.norm();
            while (scale >= Precision(1e-8)) {
                const Vec6 candidate = x + scale * delta;
                if (candidate(5) >= Precision(0)) {
                    const Vec6 candidate_residual = finite_local_residual(
                        C, Fp_old, alpha_old, candidate,
                        shear, bulk, yield_curve
                    );
                    if (candidate_residual.norm() < residual_norm) {
                        x = candidate;
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
        converged = finite_trial_from_increment(C, Fp_old, x, shear, bulk);
    }

    const Mat3 Fp_new_inv = converged.Fp.inverse();
    result.stress = sym(Fp_new_inv * converged.Se * Fp_new_inv.transpose());

    const Mat3 Cp_new = normalize_cp(converged.Fp.transpose() * converged.Fp);
    cp_to_state(result.state, Cp_new);
    result.state[6] = alpha_old + x(5);
    return result;
}

// -----------------------------------------------------------------------------
// Numerical algorithmic tangents of the complete discrete constitutive update
// -----------------------------------------------------------------------------

Mat6 tangent_small(const Vec6& strain,
                   const State& committed,
                   Precision shear,
                   Precision bulk,
                   const YieldCurve& yield_curve) {
    Mat6 C = Mat6::Zero();

    for (Index j = 0; j < 6; ++j) {
        const Precision h = derivative_step(strain(j));
        Vec6 plus = strain;
        Vec6 minus = strain;
        plus(j) += h;
        minus(j) -= h;

        const Vec6 sp = stress_voigt(update_small_strain(
            engineering_strain_tensor(plus), committed,
            shear, bulk, yield_curve
        ).stress);
        const Vec6 sm = stress_voigt(update_small_strain(
            engineering_strain_tensor(minus), committed,
            shear, bulk, yield_curve
        ).stress);

        C.col(j) = (sp - sm) / (Precision(2) * h);
    }

    return Precision(0.5) * (C + C.transpose());
}

Mat6 tangent_finite(const Vec6& strain,
                    const State& committed,
                    Precision shear,
                    Precision bulk,
                    const YieldCurve& yield_curve) {
    Mat6 C = Mat6::Zero();
    const Vec6 s0 = stress_voigt(update_finite_strain(
        engineering_strain_tensor(strain), committed,
        shear, bulk, yield_curve
    ).stress);

    for (Index j = 0; j < 6; ++j) {
        const Precision h = derivative_step(strain(j));
        Vec6 plus = strain;
        Vec6 minus = strain;
        plus(j) += h;
        minus(j) -= h;

        const bool plus_ok = admissible_green_lagrange(plus);
        const bool minus_ok = admissible_green_lagrange(minus);

        logging::error(plus_ok || minus_ok,
                       "J2: no admissible Green-Lagrange perturbation for tangent component ", j);

        if (plus_ok && minus_ok) {
            const Vec6 sp = stress_voigt(update_finite_strain(
                engineering_strain_tensor(plus), committed,
                shear, bulk, yield_curve
            ).stress);
            const Vec6 sm = stress_voigt(update_finite_strain(
                engineering_strain_tensor(minus), committed,
                shear, bulk, yield_curve
            ).stress);
            C.col(j) = (sp - sm) / (Precision(2) * h);
        } else if (plus_ok) {
            const Vec6 sp = stress_voigt(update_finite_strain(
                engineering_strain_tensor(plus), committed,
                shear, bulk, yield_curve
            ).stress);
            C.col(j) = (sp - s0) / h;
        } else {
            const Vec6 sm = stress_voigt(update_finite_strain(
                engineering_strain_tensor(minus), committed,
                shear, bulk, yield_curve
            ).stress);
            C.col(j) = (s0 - sm) / h;
        }
    }

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

template<bool Finite>
Update3D solve_shell_plane_stress(const Vec5& shell_strain,
                                  const State& committed,
                                  Precision youngs,
                                  Precision poisson,
                                  Precision shear,
                                  Precision bulk,
                                  const YieldCurve& yield_curve,
                                  Vec6& converged_strain) {
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
            current = update_finite_strain(engineering_strain_tensor(e), committed,
                                           shear, bulk, yield_curve);
        } else {
            current = update_small_strain(engineering_strain_tensor(e), committed,
                                          shear, bulk, yield_curve);
        }

        const Precision residual = current.stress(2, 2);
        if (std::abs(residual) <= tolerance) {
            converged = true;
            converged_strain = e;
            break;
        }

        const Precision h = derivative_step(e33);
        Vec6 ep = shell_to_volume_strain(shell_strain, e33 + h);
        Vec6 em = shell_to_volume_strain(shell_strain, e33 - h);

        Precision derivative = Precision(0);
        if constexpr (Finite) {
            const bool plus_ok = admissible_green_lagrange(ep);
            const bool minus_ok = admissible_green_lagrange(em);
            logging::error(plus_ok || minus_ok,
                           "J2: shell plane-stress derivative has no admissible perturbation");

            if (plus_ok && minus_ok) {
                const Precision rp = update_finite_strain(engineering_strain_tensor(ep), committed,
                                                          shear, bulk, yield_curve).stress(2, 2);
                const Precision rm = update_finite_strain(engineering_strain_tensor(em), committed,
                                                          shear, bulk, yield_curve).stress(2, 2);
                derivative = (rp - rm) / (Precision(2) * h);
            } else if (plus_ok) {
                const Precision rp = update_finite_strain(engineering_strain_tensor(ep), committed,
                                                          shear, bulk, yield_curve).stress(2, 2);
                derivative = (rp - residual) / h;
            } else {
                const Precision rm = update_finite_strain(engineering_strain_tensor(em), committed,
                                                          shear, bulk, yield_curve).stress(2, 2);
                derivative = (residual - rm) / h;
            }
        } else {
            const Precision rp = update_small_strain(engineering_strain_tensor(ep), committed,
                                                     shear, bulk, yield_curve).stress(2, 2);
            const Precision rm = update_small_strain(engineering_strain_tensor(em), committed,
                                                     shear, bulk, yield_curve).stress(2, 2);
            derivative = (rp - rm) / (Precision(2) * h);
        }

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

template<bool Finite>
Update3D solve_axial_stress(const Precision axial_strain,
                            const State& committed,
                            Precision youngs,
                            Precision poisson,
                            Precision shear,
                            Precision bulk,
                            const YieldCurve& yield_curve,
                            Vec6& converged_strain) {
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
            current = update_finite_strain(engineering_strain_tensor(e), committed,
                                           shear, bulk, yield_curve);
        } else {
            current = update_small_strain(engineering_strain_tensor(e), committed,
                                          shear, bulk, yield_curve);
        }

        Eigen::Matrix<Precision, 2, 1> residual;
        residual << current.stress(1, 1), current.stress(2, 2);

        if (residual.cwiseAbs().maxCoeff() <= tolerance) {
            converged = true;
            converged_strain = e;
            break;
        }

        Eigen::Matrix<Precision, 2, 2> J;
        for (Index column = 0; column < 2; ++column) {
            const Index component = column + 1;
            const Precision h = derivative_step(e(component));
            Vec6 ep = e;
            Vec6 em = e;
            ep(component) += h;
            em(component) -= h;

            Eigen::Matrix<Precision, 2, 1> rp;
            Eigen::Matrix<Precision, 2, 1> rm;

            if constexpr (Finite) {
                logging::error(admissible_green_lagrange(ep) && admissible_green_lagrange(em),
                               "J2: axial transverse derivative left admissible domain");
                const Update3D up = update_finite_strain(engineering_strain_tensor(ep), committed,
                                                         shear, bulk, yield_curve);
                const Update3D um = update_finite_strain(engineering_strain_tensor(em), committed,
                                                         shear, bulk, yield_curve);
                rp << up.stress(1, 1), up.stress(2, 2);
                rm << um.stress(1, 1), um.stress(2, 2);
            } else {
                const Update3D up = update_small_strain(engineering_strain_tensor(ep), committed,
                                                        shear, bulk, yield_curve);
                const Update3D um = update_small_strain(engineering_strain_tensor(em), committed,
                                                        shear, bulk, yield_curve);
                rp << up.stress(1, 1), up.stress(2, 2);
                rm << um.stress(1, 1), um.stress(2, 2);
            }

            J.col(column) = (rp - rm) / (Precision(2) * h);
        }

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
    const Update3D update = update_small_strain(strain.tensor(), committed,
                                                shear, bulk,
                                                yield_points_);

    stress = VolumeStressCauchy(update.stress);
    tangent = tangent_small(strain.voigt(), committed,
                            shear, bulk,
                            yield_points_);
    store_state(state, update.state);
}

void IsotropicJ2Elasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                     Precision*                       state,
                                     VolumeStressPK2&                 stress,
                                     Mat6&                            tangent) const {
    const State committed = load_state(state);
    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();
    const Update3D update = update_finite_strain(strain.tensor(), committed,
                                                 shear, bulk,
                                                 yield_points_);

    stress = VolumeStressPK2(update.stress);
    tangent = tangent_finite(strain.voigt(), committed,
                             shear, bulk,
                             yield_points_);
    store_state(state, update.state);
}

void IsotropicJ2Elasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                     Precision*                           state,
                                     ShellMaterialStressCauchy&            stress,
                                     Mat5&                                 tangent) const {
    const State committed = load_state(state);
    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();
    Vec6 volume_strain;
    const Update3D update = solve_shell_plane_stress<false>(
        strain.values(), committed,
        youngs, poisson, shear, bulk,
        yield_points_,
        volume_strain
    );

    const Mat6 C3 = tangent_small(volume_strain, committed,
                                  shear, bulk,
                                  yield_points_);
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
    const Update3D update = solve_shell_plane_stress<true>(
        strain.values(), committed,
        youngs, poisson, shear, bulk,
        yield_points_,
        volume_strain
    );

    const Mat6 C3 = tangent_finite(volume_strain, committed,
                                   shear, bulk,
                                   yield_points_);
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
    const Update3D update = solve_axial_stress<false>(
        strain.value(), committed,
        youngs, poisson, shear, bulk,
        yield_points_,
        volume_strain
    );

    const Mat6 C3 = tangent_small(volume_strain, committed,
                                  shear, bulk,
                                  yield_points_);
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
    const Update3D update = solve_axial_stress<true>(
        strain.value(), committed,
        youngs, poisson, shear, bulk,
        yield_points_,
        volume_strain
    );

    const Mat6 C3 = tangent_finite(volume_strain, committed,
                                   shear, bulk,
                                   yield_points_);
    stress.value() = update.stress(0, 0);
    tangent = condense_axial_tangent(C3);
    store_state(state, update.state);
}

} // namespace fem::material
