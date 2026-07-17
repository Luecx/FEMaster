#include "neo_hooke_elasticity.h"

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

#include <Eigen/LU>

#include <array>
#include <cmath>

namespace fem::material {

NeoHookeElasticity::NeoHookeElasticity(Precision youngs_in, Precision poisson_in)
    : youngs     (youngs_in),
      poisson    (poisson_in),
      lame_lambda(youngs_in * poisson_in /
                  ((Precision(1) + poisson_in) * (Precision(1) - Precision(2) * poisson_in))),
      mu         (youngs_in / (Precision(2) * (Precision(1) + poisson_in))) {
    logging::error(youngs > Precision(0),
                   "NEO_HOOKE: Young's modulus must be positive");
    logging::error(poisson > Precision(-1) && poisson < Precision(0.5),
                   "NEO_HOOKE: Poisson ratio must be in (-1, 0.5)");
}

bool NeoHookeElasticity::supports_axial_linearized() const {
    return true;
}

bool NeoHookeElasticity::supports_axial_green_lagrange() const {
    return true;
}

bool NeoHookeElasticity::supports_volume_linearized() const {
    return true;
}

bool NeoHookeElasticity::supports_volume_green_lagrange() const {
    return true;
}

bool NeoHookeElasticity::supports_shell_integration_linearized() const {
    return true;
}

bool NeoHookeElasticity::supports_shell_integration_green_lagrange() const {
    return true;
}

Mat6 NeoHookeElasticity::linear_tangent() const {
    Mat6 tangent;
    Precision C11 = lame_lambda + Precision(2) * mu;
    tangent <<
        C11         , lame_lambda , lame_lambda , Precision(0), Precision(0), Precision(0),
        lame_lambda , C11         , lame_lambda , Precision(0), Precision(0), Precision(0),
        lame_lambda , lame_lambda , C11         , Precision(0), Precision(0), Precision(0),
        Precision(0), Precision(0), Precision(0), mu,           Precision(0), Precision(0),
        Precision(0), Precision(0), Precision(0), Precision(0), mu,           Precision(0),
        Precision(0), Precision(0), Precision(0), Precision(0), Precision(0), mu;
    return tangent;
}

Mat5 NeoHookeElasticity::linear_shell_tangent() const {
    const Precision scalar = youngs / (Precision(1) - poisson * poisson);

    Mat5 tangent = Mat5::Zero();
    tangent(0, 0) = scalar;
    tangent(0, 1) = scalar * poisson;
    tangent(1, 0) = scalar * poisson;
    tangent(1, 1) = scalar;
    tangent(2, 2) = mu;
    tangent(3, 3) = mu;
    tangent(4, 4) = mu;
    return tangent;
}

void NeoHookeElasticity::evaluate(const AxialStrainLinearized& strain,
                                  const Precision*             state_old,
                                  Precision*                   state_new,
                                  AxialStressCauchy&           stress,
                                  Precision&                   tangent) const {
    (void) state_old;
    (void) state_new;

    tangent        = youngs;
    stress.value() = tangent * strain.value();
}

void NeoHookeElasticity::evaluate(const AxialStrainGreenLagrange& strain,
                                  const Precision*                state_old,
                                  Precision*                      state_new,
                                  AxialStressPK2&                 stress,
                                  Precision&                      tangent) const {
    (void) state_old;
    (void) state_new;

    const Precision c = Precision(1) + Precision(2) * strain.value();

    logging::error(c > Precision(0),
                   "NEO_HOOKE: non-positive axial right Cauchy-Green component");

    Precision y = -poisson * std::log(c);

    for (Index i = 0; i < 20; ++i) {
        const Precision x  = std::exp(y);
        const Precision r  = mu * (x - Precision(1))
                           + lame_lambda * (Precision(0.5) * std::log(c) + y);
        const Precision dy = -r / (mu * x + lame_lambda);

        y += dy;

        if (std::abs(dy) < Precision(1e-12)) {
            break;
        }
    }

    const Precision x = std::exp(y);

    stress.value() = mu * (Precision(1) - x / c);
    tangent        = mu * x / (c * c)
                   * (Precision(2) * mu * x + Precision(3) * lame_lambda)
                   / (mu * x + lame_lambda);
}

void NeoHookeElasticity::evaluate(const VolumeStrainLinearized& strain,
                                  const Precision*              state_old,
                                  Precision*                    state_new,
                                  VolumeStressCauchy&           stress,
                                  Mat6&                         tangent) const {
    (void) state_old;
    (void) state_new;

    tangent        = linear_tangent();
    stress.voigt() = tangent * strain.voigt();
}

void NeoHookeElasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                  const Precision*                 state_old,
                                  Precision*                       state_new,
                                  VolumeStressPK2&                 stress,
                                  Mat6&                            tangent) const {
    (void) state_old;
    (void) state_new;

    const Mat3 E    = strain.tensor();
    const Mat3 C    = Mat3::Identity() + Precision(2) * E;
    const Precision det_c = C.determinant();

    logging::error(det_c > Precision(0),
                   "NEO_HOOKE: non-positive right Cauchy-Green determinant");

    const Mat3      C_inv = C.inverse();
    const Precision log_j = Precision(0.5) * std::log(det_c);

    const Mat3 S =
        mu * (Mat3::Identity() - C_inv)
      + lame_lambda * log_j * C_inv;

    stress = VolumeStressPK2(S);

    const Precision q = mu - lame_lambda * log_j;

    tangent.setZero();
    for (Index col = 0; col < 6; ++col) {
        Mat3 dE = Mat3::Zero();

        if (col < 3) {
            dE(col, col) = Precision(1);
        } else if (col == 3) {
            dE(1, 2) = Precision(0.5);
            dE(2, 1) = Precision(0.5);
        } else if (col == 4) {
            dE(0, 2) = Precision(0.5);
            dE(2, 0) = Precision(0.5);
        } else {
            dE(0, 1) = Precision(0.5);
            dE(1, 0) = Precision(0.5);
        }

        const Mat3 dS =
            lame_lambda * (C_inv * dE).trace() * C_inv
          + Precision(2) * q * C_inv * dE * C_inv;

        tangent.col(col) = VolumeStressPK2(dS).voigt();
    }
}

void NeoHookeElasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                  const Precision*                     state_old,
                                  Precision*                           state_new,
                                  ShellMaterialStressCauchy&            stress,
                                  Mat5&                                 tangent) const {
    (void) state_old;
    (void) state_new;

    tangent         = linear_shell_tangent();
    stress.values() = tangent * strain.values();
}

void NeoHookeElasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                                  const Precision*                        state_old,
                                  Precision*                              state_new,
                                  ShellMaterialStressPK2&                 stress,
                                  Mat5&                                   tangent) const {
    (void) state_old;
    (void) state_new;

    Mat3 C = Mat3::Identity();
    C(0, 0) = Precision(1) + Precision(2) * strain.values()(0);
    C(1, 1) = Precision(1) + Precision(2) * strain.values()(1);
    C(0, 1) = strain.values()(2);
    C(1, 0) = strain.values()(2);
    C(0, 2) = strain.values()(3);
    C(2, 0) = strain.values()(3);
    C(1, 2) = strain.values()(4);
    C(2, 1) = strain.values()(4);

    const Mat2 in_plane = C.template block<2, 2>(0, 0);
    const Vec2 shear_column(C(0, 2), C(1, 2));
    const Precision det_in_plane = in_plane.determinant();

    logging::error(det_in_plane > Precision(0),
                   "NEO_HOOKE: non-positive shell in-plane right Cauchy-Green determinant");

    const Precision min_c33 = (shear_column.transpose() * in_plane.inverse() * shear_column)(0, 0);
    Precision c33 = Precision(1) - Precision(2) * poisson * (strain.values()(0) + strain.values()(1));

    if (c33 <= min_c33) {
        c33 = min_c33 + Precision(1e-10);
    }

    auto evaluate_full = [&](const Mat3& c, Mat3& full_stress, Mat6& full_tangent) {
        const Precision det_c = c.determinant();

        logging::error(det_c > Precision(0),
                       "NEO_HOOKE: non-positive shell right Cauchy-Green determinant");

        const Mat3      c_inv = c.inverse();
        const Precision log_j = Precision(0.5) * std::log(det_c);

        full_stress =
            mu * (Mat3::Identity() - c_inv)
          + lame_lambda * log_j * c_inv;

        const Precision q = mu - lame_lambda * log_j;

        full_tangent.setZero();
        for (Index col = 0; col < 6; ++col) {
            Mat3 dE = Mat3::Zero();

            if (col < 3) {
                dE(col, col) = Precision(1);
            } else if (col == 3) {
                dE(1, 2) = Precision(0.5);
                dE(2, 1) = Precision(0.5);
            } else if (col == 4) {
                dE(0, 2) = Precision(0.5);
                dE(2, 0) = Precision(0.5);
            } else {
                dE(0, 1) = Precision(0.5);
                dE(1, 0) = Precision(0.5);
            }

            const Mat3 dS =
                lame_lambda * (c_inv * dE).trace() * c_inv
              + Precision(2) * q * c_inv * dE * c_inv;

            full_tangent.col(col) = VolumeStressPK2(dS).voigt();
        }
    };

    Mat3 full_stress;
    Mat6 full_tangent;

    for (Index i = 0; i < 30; ++i) {
        C(2, 2) = c33;
        evaluate_full(C, full_stress, full_tangent);

        const Precision residual = full_stress(2, 2);
        const Precision slope    = Precision(0.5) * full_tangent(2, 2);

        logging::error(std::abs(slope) > Precision(1e-20),
                       "NEO_HOOKE: singular shell plane-stress thickness iteration");

        Precision dc33 = -residual / slope;
        c33 += dc33;

        if (c33 <= min_c33) {
            c33 = Precision(0.5) * (C(2, 2) + min_c33 + Precision(1e-10));
            dc33 = c33 - C(2, 2);
        }

        if (std::abs(dc33) < Precision(1e-12)) {
            break;
        }
    }

    C(2, 2) = c33;
    evaluate_full(C, full_stress, full_tangent);

    stress.values() << full_stress(0, 0),
                       full_stress(1, 1),
                       full_stress(0, 1),
                       full_stress(0, 2),
                       full_stress(1, 2);

    const std::array<Index, 5> shell_rows { 0, 1, 5, 4, 3 };
    const std::array<Index, 5> shell_cols { 0, 1, 5, 4, 3 };

    tangent.setZero();
    for (Index row = 0; row < 5; ++row) {
        for (Index col = 0; col < 5; ++col) {
            tangent(row, col) =
                full_tangent(shell_rows[row], shell_cols[col])
              - full_tangent(shell_rows[row], 2)
              * full_tangent(2, shell_cols[col])
              / full_tangent(2, 2);
        }
    }
}

} // namespace fem::material
