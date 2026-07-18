#include "isotropic_elasticity.h"

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

namespace fem::material {

IsotropicElasticity::IsotropicElasticity(Precision youngs_in, Precision poisson_in)
    : youngs (youngs_in),
      poisson(poisson_in),
      shear  (youngs_in / (Precision(2) * (Precision(1) + poisson_in))) {
    logging::error(youngs > Precision(0),
                   "ISOTROPIC: Young's modulus must be positive");
    logging::error(poisson > Precision(-1) && poisson < Precision(0.5),
                   "ISOTROPIC: Poisson ratio must be in (-1, 0.5)");
}

bool IsotropicElasticity::supports_axial_linearized() const {
    return true;
}

bool IsotropicElasticity::supports_axial_green_lagrange() const {
    return true;
}

bool IsotropicElasticity::supports_volume_linearized() const {
    return true;
}

bool IsotropicElasticity::supports_volume_green_lagrange() const {
    return true;
}

bool IsotropicElasticity::supports_shell_integration_linearized() const {
    return true;
}

bool IsotropicElasticity::supports_shell_integration_green_lagrange() const {
    return true;
}

Mat3 IsotropicElasticity::plane_stress_tangent() const {
    const Precision scalar = youngs / (Precision(1) - poisson * poisson);
    Mat3 tangent;
    tangent << scalar,           scalar * poisson, Precision(0),
               scalar * poisson, scalar,           Precision(0),
               Precision(0),     Precision(0),     shear;
    return tangent;
}

Mat5 IsotropicElasticity::shell_material_tangent() const {
    Mat5 tangent = Mat5::Zero();
    tangent.template block<3, 3>(0, 0) = plane_stress_tangent();
    tangent(3, 3) = shear;
    tangent(4, 4) = shear;
    return tangent;
}

Mat6 IsotropicElasticity::volume_tangent() const {
    const Precision scalar = youngs
        / ((Precision(1) + poisson) * (Precision(1) - Precision(2) * poisson));
    const Precision mu = Precision(1) - Precision(2) * poisson;

    Mat6 tangent;
    tangent <<
        Precision(1) - poisson, poisson, poisson, Precision(0), Precision(0), Precision(0),
        poisson, Precision(1) - poisson, poisson, Precision(0), Precision(0), Precision(0),
        poisson, poisson, Precision(1) - poisson, Precision(0), Precision(0), Precision(0),
        Precision(0), Precision(0), Precision(0), mu / Precision(2), Precision(0), Precision(0),
        Precision(0), Precision(0), Precision(0), Precision(0), mu / Precision(2), Precision(0),
        Precision(0), Precision(0), Precision(0), Precision(0), Precision(0), mu / Precision(2);
    return scalar * tangent;
}

void IsotropicElasticity::evaluate(const AxialStrainLinearized& strain,
                                   const Precision*             state_old,
                                   Precision*                   state_new,
                                   AxialStressCauchy&           stress,
                                   Precision&                   tangent) const {
    (void) state_old;
    (void) state_new;
    tangent        = youngs;
    stress.value() = tangent * strain.value();
}

void IsotropicElasticity::evaluate(const AxialStrainGreenLagrange& strain,
                                   const Precision*                state_old,
                                   Precision*                      state_new,
                                   AxialStressPK2&                 stress,
                                   Precision&                      tangent) const {
    (void) state_old;
    (void) state_new;
    tangent        = youngs;
    stress.value() = tangent * strain.value();
}

void IsotropicElasticity::evaluate(const VolumeStrainLinearized& strain,
                                   const Precision*              state_old,
                                   Precision*                    state_new,
                                   VolumeStressCauchy&           stress,
                                   Mat6&                         tangent) const {
    (void) state_old;
    (void) state_new;
    tangent        = volume_tangent();
    stress.voigt() = tangent * strain.voigt();
}

void IsotropicElasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                   const Precision*                 state_old,
                                   Precision*                       state_new,
                                   VolumeStressPK2&                 stress,
                                   Mat6&                            tangent) const {
    (void) state_old;
    (void) state_new;
    tangent        = volume_tangent();
    stress.voigt() = tangent * strain.voigt();
}

void IsotropicElasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                   const Precision*                     state_old,
                                   Precision*                           state_new,
                                   ShellMaterialStressCauchy&            stress,
                                   Mat5&                                 tangent) const {
    (void) state_old;
    (void) state_new;
    tangent         = shell_material_tangent();
    stress.values() = tangent * strain.values();
}

void IsotropicElasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                                   const Precision*                        state_old,
                                   Precision*                              state_new,
                                   ShellMaterialStressPK2&                 stress,
                                   Mat5&                                   tangent) const {
    (void) state_old;
    (void) state_new;
    tangent         = shell_material_tangent();
    stress.values() = tangent * strain.values();
}

} // namespace fem::material
