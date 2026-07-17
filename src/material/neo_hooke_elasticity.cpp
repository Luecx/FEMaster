#include "neo_hooke_elasticity.h"

#include "../core/logging.h"
#include "strain/axial_strain_green_lagrange.h"
#include "strain/axial_strain_linearized.h"
#include "strain/volume_strain_green_lagrange.h"
#include "strain/volume_strain_linearized.h"
#include "stress/axial_stress_cauchy.h"
#include "stress/axial_stress_pk2.h"
#include "stress/volume_stress_cauchy.h"
#include "stress/volume_stress_pk2.h"

#include <Eigen/LU>

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

} // namespace fem::material
