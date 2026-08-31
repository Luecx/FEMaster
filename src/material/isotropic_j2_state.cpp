/**
 * @file isotropic_j2_state.cpp
 * @brief Implements immutable source/target state semantics for J2 plasticity.
 *
 * The plastic return-mapping kernels in isotropic_j2_elasticity.cpp operate on a
 * mutable working row. This translation unit confines that mutation to the
 * target state supplied to `integrate()`. The public source state is never
 * modified.
 *
 * Read-only `evaluate()` is deliberately different from `integrate()`: it
 * evaluates stress and tangent with the supplied plastic history frozen. This is
 * the required operation for result recovery and diagnostics because evaluating
 * an already accepted configuration must not perform another plastic increment.
 *
 * @see IsotropicJ2Elasticity
 * @see Elasticity
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

namespace fem::material {
namespace {

constexpr Index j2_state_size = 7;

using State = std::array<Precision, j2_state_size>;

State load_state_readonly(const Precision* state) {
    logging::error(state != nullptr, "J2: source material state is null");

    State values{};
    for (Index i = 0; i < j2_state_size; ++i) {
        values[static_cast<std::size_t>(i)] = state[i];
    }
    return values;
}

Mat3 sym(const Mat3& A) {
    return Precision(0.5) * (A + A.transpose());
}

Mat3 dev(const Mat3& A) {
    return A - A.trace() / Precision(3) * Mat3::Identity();
}

Mat3 cp_from_state(const State& state) {
    Mat3 Cp;
    Cp << state[0], state[5], state[4],
          state[5], state[1], state[3],
          state[4], state[3], state[2];
    return sym(Cp);
}

template<typename Function>
Mat3 spd_function(const Mat3& A, Function&& function, const char* name) {
    Eigen::SelfAdjointEigenSolver<Mat3> solver(sym(A));
    logging::error(solver.info() == Eigen::Success,
                   "J2: eigendecomposition failed for ", name);

    const Vec3 eigenvalues = solver.eigenvalues();
    const Precision scale = std::max(Precision(1), eigenvalues.cwiseAbs().maxCoeff());
    logging::error(
        eigenvalues.minCoeff() > Precision(100) * std::numeric_limits<Precision>::epsilon() * scale,
        "J2: ", name, " is not symmetric positive definite"
    );

    Vec3 mapped;
    for (Index i = 0; i < 3; ++i) {
        mapped(i) = function(eigenvalues(i));
    }

    return sym(solver.eigenvectors() * mapped.asDiagonal() * solver.eigenvectors().transpose());
}

Mat3 spd_sqrt(const Mat3& A) {
    return spd_function(A, [](Precision x) { return std::sqrt(x); }, "plastic metric");
}

Mat3 spd_log(const Mat3& A) {
    return spd_function(A, [](Precision x) { return std::log(x); }, "plastic metric");
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

Precision derivative_step(Precision value) {
    return Precision(1e-7) * std::max(Precision(1), std::abs(value));
}

Mat6 elastic_volume_tangent(Precision shear, Precision bulk) {
    const Precision lame = bulk - Precision(2) * shear / Precision(3);
    Mat6 C = Mat6::Zero();

    for (Index i = 0; i < 3; ++i) {
        for (Index j = 0; j < 3; ++j) {
            C(i, j) = lame;
        }
        C(i, i) += Precision(2) * shear;
    }

    C(3, 3) = shear;
    C(4, 4) = shear;
    C(5, 5) = shear;
    return C;
}

Mat3 fixed_small_stress(const Vec6& strain,
                        const State& state,
                        Precision shear,
                        Precision bulk) {
    const Mat3 plastic_strain = Precision(0.5) * spd_log(cp_from_state(state));
    const Mat3 elastic_strain = sym(engineering_strain_tensor(strain) - plastic_strain);

    return sym(
        bulk * elastic_strain.trace() * Mat3::Identity()
      + Precision(2) * shear * dev(elastic_strain)
    );
}

bool admissible_green_lagrange(const Vec6& strain) {
    const Mat3 C = Mat3::Identity() + Precision(2) * engineering_strain_tensor(strain);
    Eigen::SelfAdjointEigenSolver<Mat3> solver(sym(C), Eigen::EigenvaluesOnly);
    return solver.info() == Eigen::Success
        && solver.eigenvalues().minCoeff() > Precision(1e-12);
}

Mat3 fixed_finite_stress(const Vec6& strain,
                         const State& state,
                         Precision shear,
                         Precision bulk) {
    const Mat3 C = sym(Mat3::Identity() + Precision(2) * engineering_strain_tensor(strain));
    logging::error(admissible_green_lagrange(strain),
                   "J2: read-only finite evaluation received inadmissible Green-Lagrange strain");

    const Mat3 Fp = spd_sqrt(cp_from_state(state));
    const Mat3 Fp_inv = Fp.inverse();
    const Mat3 Ce = sym(Fp_inv.transpose() * C * Fp_inv);

    // Validate the elastic metric before constructing the StVK response.
    (void) spd_log(Ce);

    const Mat3 Ee = Precision(0.5) * (Ce - Mat3::Identity());
    const Precision lame = bulk - Precision(2) * shear / Precision(3);
    const Mat3 Se = lame * Ee.trace() * Mat3::Identity()
                  + Precision(2) * shear * Ee;

    return sym(Fp_inv * Se * Fp_inv.transpose());
}

Mat6 fixed_finite_tangent(const Vec6& strain,
                          const State& state,
                          Precision shear,
                          Precision bulk) {
    Mat6 C = Mat6::Zero();
    const Vec6 s0 = stress_voigt(fixed_finite_stress(strain, state, shear, bulk));

    for (Index j = 0; j < 6; ++j) {
        const Precision h = derivative_step(strain(j));
        Vec6 plus  = strain;
        Vec6 minus = strain;
        plus(j)  += h;
        minus(j) -= h;

        const bool plus_ok  = admissible_green_lagrange(plus);
        const bool minus_ok = admissible_green_lagrange(minus);
        logging::error(plus_ok || minus_ok,
                       "J2: no admissible frozen-state tangent perturbation for component ", j);

        if (plus_ok && minus_ok) {
            const Vec6 sp = stress_voigt(fixed_finite_stress(plus, state, shear, bulk));
            const Vec6 sm = stress_voigt(fixed_finite_stress(minus, state, shear, bulk));
            C.col(j) = (sp - sm) / (Precision(2) * h);
        } else if (plus_ok) {
            const Vec6 sp = stress_voigt(fixed_finite_stress(plus, state, shear, bulk));
            C.col(j) = (sp - s0) / h;
        } else {
            const Vec6 sm = stress_voigt(fixed_finite_stress(minus, state, shear, bulk));
            C.col(j) = (s0 - sm) / h;
        }
    }

    return Precision(0.5) * (C + C.transpose());
}

template<bool Finite>
Mat3 fixed_volume_stress(const Vec6& strain,
                         const State& state,
                         Precision shear,
                         Precision bulk) {
    if constexpr (Finite) {
        return fixed_finite_stress(strain, state, shear, bulk);
    } else {
        return fixed_small_stress(strain, state, shear, bulk);
    }
}

template<bool Finite>
Vec6 solve_axial_frozen(Precision axial_strain,
                        const State& state,
                        Precision poisson,
                        Precision shear,
                        Precision bulk) {
    Vec6 strain = Vec6::Zero();
    strain(0) = axial_strain;
    strain(1) = -poisson * axial_strain;
    strain(2) = -poisson * axial_strain;

    if constexpr (Finite) {
        strain(1) = std::max(strain(1), Precision(-0.49));
        strain(2) = std::max(strain(2), Precision(-0.49));
    }

    const Precision tolerance = Precision(1e-10)
        * std::max(Precision(1), Precision(3) * bulk);

    for (Index iteration = 0; iteration < 30; ++iteration) {
        const Mat3 stress = fixed_volume_stress<Finite>(strain, state, shear, bulk);
        Eigen::Matrix<Precision, 2, 1> residual;
        residual << stress(1, 1), stress(2, 2);

        if (residual.cwiseAbs().maxCoeff() <= tolerance) {
            return strain;
        }

        Eigen::Matrix<Precision, 2, 2> J;
        for (Index column = 0; column < 2; ++column) {
            const Index component = column + 1;
            const Precision h = derivative_step(strain(component));
            Vec6 plus  = strain;
            Vec6 minus = strain;
            plus(component)  += h;
            minus(component) -= h;

            const Mat3 sp = fixed_volume_stress<Finite>(plus, state, shear, bulk);
            const Mat3 sm = fixed_volume_stress<Finite>(minus, state, shear, bulk);
            J(0, column) = (sp(1, 1) - sm(1, 1)) / (Precision(2) * h);
            J(1, column) = (sp(2, 2) - sm(2, 2)) / (Precision(2) * h);
        }

        Eigen::FullPivLU<Eigen::Matrix<Precision, 2, 2>> lu(J);
        logging::error(lu.isInvertible(),
                       "J2: singular frozen-state transverse Jacobian");
        const Eigen::Matrix<Precision, 2, 1> delta = lu.solve(-residual);

        Precision scale = Precision(1);
        bool accepted = false;
        while (scale >= Precision(1e-8)) {
            Vec6 candidate = strain;
            candidate(1) += scale * delta(0);
            candidate(2) += scale * delta(1);

            if constexpr (Finite) {
                if (!admissible_green_lagrange(candidate)) {
                    scale *= Precision(0.5);
                    continue;
                }
            }

            strain = candidate;
            accepted = true;
            break;
        }

        logging::error(accepted,
                       "J2: frozen-state axial reduction could not find admissible step");
    }

    logging::error(false, "J2: frozen-state axial reduction did not converge");
    return strain;
}

template<bool Finite>
Precision frozen_axial_stress(Precision axial_strain,
                              const State& state,
                              Precision poisson,
                              Precision shear,
                              Precision bulk) {
    const Vec6 strain = solve_axial_frozen<Finite>(
        axial_strain, state, poisson, shear, bulk
    );
    return fixed_volume_stress<Finite>(strain, state, shear, bulk)(0, 0);
}

template<bool Finite>
Precision frozen_axial_tangent(Precision axial_strain,
                               const State& state,
                               Precision poisson,
                               Precision shear,
                               Precision bulk) {
    const Precision h = derivative_step(axial_strain);
    const Precision sp = frozen_axial_stress<Finite>(
        axial_strain + h, state, poisson, shear, bulk
    );
    const Precision sm = frozen_axial_stress<Finite>(
        axial_strain - h, state, poisson, shear, bulk
    );
    return (sp - sm) / (Precision(2) * h);
}

Vec6 shell_to_volume(const Vec5& shell, Precision e33) {
    Vec6 strain;
    strain << shell(0), shell(1), e33, shell(4), shell(3), shell(2);
    return strain;
}

Vec5 volume_to_shell_stress(const Mat3& stress) {
    Vec5 result;
    result << stress(0, 0), stress(1, 1), stress(0, 1), stress(0, 2), stress(1, 2);
    return result;
}

template<bool Finite>
Precision solve_shell_e33_frozen(const Vec5& shell,
                                 const State& state,
                                 Precision poisson,
                                 Precision shear,
                                 Precision bulk) {
    Precision e33 = -poisson / (Precision(1) - poisson) * (shell(0) + shell(1));
    if constexpr (Finite) {
        e33 = std::max(e33, Precision(-0.49));
    }

    const Precision tolerance = Precision(1e-10)
        * std::max(Precision(1), Precision(3) * bulk);

    for (Index iteration = 0; iteration < 30; ++iteration) {
        const Vec6 strain = shell_to_volume(shell, e33);
        const Mat3 stress = fixed_volume_stress<Finite>(strain, state, shear, bulk);
        const Precision residual = stress(2, 2);

        if (std::abs(residual) <= tolerance) {
            return e33;
        }

        const Precision h = derivative_step(e33);
        const Vec6 plus  = shell_to_volume(shell, e33 + h);
        const Vec6 minus = shell_to_volume(shell, e33 - h);

        Precision derivative = Precision(0);
        if constexpr (Finite) {
            const bool plus_ok  = admissible_green_lagrange(plus);
            const bool minus_ok = admissible_green_lagrange(minus);
            logging::error(plus_ok || minus_ok,
                           "J2: shell frozen-state derivative has no admissible perturbation");

            if (plus_ok && minus_ok) {
                derivative = (
                    fixed_finite_stress(plus, state, shear, bulk)(2, 2)
                  - fixed_finite_stress(minus, state, shear, bulk)(2, 2)
                ) / (Precision(2) * h);
            } else if (plus_ok) {
                derivative = (
                    fixed_finite_stress(plus, state, shear, bulk)(2, 2) - residual
                ) / h;
            } else {
                derivative = (
                    residual - fixed_finite_stress(minus, state, shear, bulk)(2, 2)
                ) / h;
            }
        } else {
            derivative = (
                fixed_small_stress(plus, state, shear, bulk)(2, 2)
              - fixed_small_stress(minus, state, shear, bulk)(2, 2)
            ) / (Precision(2) * h);
        }

        logging::error(std::isfinite(derivative)
                       && std::abs(derivative) > Precision(100) * std::numeric_limits<Precision>::epsilon(),
                       "J2: singular frozen-state shell thickness derivative");

        const Precision delta = -residual / derivative;
        Precision scale = Precision(1);
        bool accepted = false;

        while (scale >= Precision(1e-8)) {
            const Precision candidate_e33 = e33 + scale * delta;
            const Vec6 candidate = shell_to_volume(shell, candidate_e33);

            if constexpr (Finite) {
                if (!admissible_green_lagrange(candidate)) {
                    scale *= Precision(0.5);
                    continue;
                }
            }

            e33 = candidate_e33;
            accepted = true;
            break;
        }

        logging::error(accepted,
                       "J2: frozen-state shell reduction could not find admissible step");
    }

    logging::error(false, "J2: frozen-state shell reduction did not converge");
    return e33;
}

template<bool Finite>
Vec5 frozen_shell_stress(const Vec5& shell,
                         const State& state,
                         Precision poisson,
                         Precision shear,
                         Precision bulk) {
    const Precision e33 = solve_shell_e33_frozen<Finite>(
        shell, state, poisson, shear, bulk
    );
    return volume_to_shell_stress(
        fixed_volume_stress<Finite>(shell_to_volume(shell, e33), state, shear, bulk)
    );
}

template<bool Finite>
Mat5 frozen_shell_tangent(const Vec5& shell,
                          const State& state,
                          Precision poisson,
                          Precision shear,
                          Precision bulk) {
    Mat5 tangent = Mat5::Zero();

    for (Index j = 0; j < 5; ++j) {
        const Precision h = derivative_step(shell(j));
        Vec5 plus  = shell;
        Vec5 minus = shell;
        plus(j)  += h;
        minus(j) -= h;

        const Vec5 sp = frozen_shell_stress<Finite>(plus, state, poisson, shear, bulk);
        const Vec5 sm = frozen_shell_stress<Finite>(minus, state, poisson, shear, bulk);
        tangent.col(j) = (sp - sm) / (Precision(2) * h);
    }

    return Precision(0.5) * (tangent + tangent.transpose());
}

void copy_source_to_target(const Precision* state, Precision* target_state) {
    logging::error(state != nullptr, "J2: source material state is null");
    logging::error(target_state != nullptr, "J2: target material state is null");
    logging::error(state != target_state,
                   "J2: source and target material state must not alias");

    for (Index i = 0; i < j2_state_size; ++i) {
        target_state[i] = state[i];
    }
}

} // namespace

// -----------------------------------------------------------------------------
// Read-only response
// -----------------------------------------------------------------------------

void IsotropicJ2Elasticity::evaluate(const VolumeStrainLinearized& strain,
                                     const Precision*              state,
                                     VolumeStressCauchy&           stress,
                                     Mat6&                         tangent) const {
    const State history = load_state_readonly(state);
    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();

    stress = VolumeStressCauchy(fixed_small_stress(strain.voigt(), history, shear, bulk));
    tangent = elastic_volume_tangent(shear, bulk);
}

void IsotropicJ2Elasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                     const Precision*                   state,
                                     VolumeStressPK2&                 stress,
                                     Mat6&                           tangent) const {
    const State history = load_state_readonly(state);
    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();

    stress = VolumeStressPK2(fixed_finite_stress(strain.voigt(), history, shear, bulk));
    tangent = fixed_finite_tangent(strain.voigt(), history, shear, bulk);
}

void IsotropicJ2Elasticity::evaluate(const AxialStrainLinearized& strain,
                                     const Precision*             state,
                                     AxialStressCauchy&           stress,
                                     Precision&                   tangent) const {
    const State history = load_state_readonly(state);
    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();

    stress.value() = frozen_axial_stress<false>(
        strain.value(), history, poisson, shear, bulk
    );
    tangent = frozen_axial_tangent<false>(
        strain.value(), history, poisson, shear, bulk
    );
}

void IsotropicJ2Elasticity::evaluate(const AxialStrainGreenLagrange& strain,
                                     const Precision*                  state,
                                     AxialStressPK2&                 stress,
                                     Precision&                     tangent) const {
    const State history = load_state_readonly(state);
    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();

    stress.value() = frozen_axial_stress<true>(
        strain.value(), history, poisson, shear, bulk
    );
    tangent = frozen_axial_tangent<true>(
        strain.value(), history, poisson, shear, bulk
    );
}

void IsotropicJ2Elasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                     const Precision*                     state,
                                     ShellMaterialStressCauchy&            stress,
                                     Mat5&                               tangent) const {
    const State history = load_state_readonly(state);
    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();

    stress.values() = frozen_shell_stress<false>(
        strain.values(), history, poisson, shear, bulk
    );
    tangent = frozen_shell_tangent<false>(
        strain.values(), history, poisson, shear, bulk
    );
}

void IsotropicJ2Elasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                                     const Precision*                          state,
                                     ShellMaterialStressPK2&                 stress,
                                     Mat5&                                  tangent) const {
    const State history = load_state_readonly(state);
    const Precision shear = shear_modulus();
    const Precision bulk  = bulk_modulus();

    stress.values() = frozen_shell_stress<true>(
        strain.values(), history, poisson, shear, bulk
    );
    tangent = frozen_shell_tangent<true>(
        strain.values(), history, poisson, shear, bulk
    );
}

// -----------------------------------------------------------------------------
// Source -> target state integration
// -----------------------------------------------------------------------------

void IsotropicJ2Elasticity::integrate(const AxialStrainLinearized& strain,
                                      const Precision*             state,
                                      Precision*                   target_state,
                                      AxialStressCauchy&           stress,
                                      Precision&                   tangent) const {
    copy_source_to_target(state, target_state);
    evaluate(strain, target_state, stress, tangent);
}

void IsotropicJ2Elasticity::integrate(const AxialStrainGreenLagrange& strain,
                                      const Precision*                  state,
                                      Precision*                      target_state,
                                      AxialStressPK2&                 stress,
                                      Precision&                     tangent) const {
    copy_source_to_target(state, target_state);
    evaluate(strain, target_state, stress, tangent);
}

void IsotropicJ2Elasticity::integrate(const VolumeStrainLinearized& strain,
                                      const Precision*              state,
                                      Precision*                    target_state,
                                      VolumeStressCauchy&           stress,
                                      Mat6&                         tangent) const {
    copy_source_to_target(state, target_state);
    evaluate(strain, target_state, stress, tangent);
}

void IsotropicJ2Elasticity::integrate(const VolumeStrainGreenLagrange& strain,
                                      const Precision*                   state,
                                      Precision*                       target_state,
                                      VolumeStressPK2&                 stress,
                                      Mat6&                           tangent) const {
    copy_source_to_target(state, target_state);
    evaluate(strain, target_state, stress, tangent);
}

void IsotropicJ2Elasticity::integrate(const ShellMaterialStrainLinearized& strain,
                                      const Precision*                     state,
                                      Precision*                           target_state,
                                      ShellMaterialStressCauchy&            stress,
                                      Mat5&                               tangent) const {
    copy_source_to_target(state, target_state);
    evaluate(strain, target_state, stress, tangent);
}

void IsotropicJ2Elasticity::integrate(const ShellMaterialStrainGreenLagrange& strain,
                                      const Precision*                          state,
                                      Precision*                              target_state,
                                      ShellMaterialStressPK2&                 stress,
                                      Mat5&                                  tangent) const {
    copy_source_to_target(state, target_state);
    evaluate(strain, target_state, stress, tangent);
}

} // namespace fem::material
