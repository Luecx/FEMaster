#include "generalised_isotropic_elasticity.h"

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

GeneralisedIsotropicElasticity::GeneralisedIsotropicElasticity(Precision youngs_in,
                                                               Precision poisson_in,
                                                               Precision shear_in)
    : youngs (youngs_in),
      poisson(poisson_in),
      shear  (shear_in) {}

bool GeneralisedIsotropicElasticity::supports_axial_linearized() const {
    return true;
}

bool GeneralisedIsotropicElasticity::supports_axial_green_lagrange() const {
    return true;
}

bool GeneralisedIsotropicElasticity::supports_volume_linearized() const {
    return true;
}

bool GeneralisedIsotropicElasticity::supports_volume_green_lagrange() const {
    return true;
}

bool GeneralisedIsotropicElasticity::supports_shell_integration_linearized() const {
    return true;
}

bool GeneralisedIsotropicElasticity::supports_shell_integration_green_lagrange() const {
    return true;
}

Mat3 GeneralisedIsotropicElasticity::plane_stress_tangent() const {
    const Precision scalar = youngs / (Precision(1) - poisson * poisson);
    Mat3 tangent;
    tangent << scalar,           scalar * poisson, Precision(0),
               scalar * poisson, scalar,           Precision(0),
               Precision(0),     Precision(0),     shear;
    return tangent;
}

Mat5 GeneralisedIsotropicElasticity::shell_material_tangent() const {
    Mat5 tangent = Mat5::Zero();
    tangent.template block<3, 3>(0, 0) = plane_stress_tangent();
    tangent(3, 3) = shear;
    tangent(4, 4) = shear;
    return tangent;
}

Mat6 GeneralisedIsotropicElasticity::volume_tangent() const {
    const Precision scalar = youngs
        / ((Precision(1) + poisson) * (Precision(1) - Precision(2) * poisson));
    const Precision c11 = (Precision(1) - poisson) * scalar;
    const Precision c12 = poisson * scalar;

    Mat6 tangent;
    tangent <<
        c11, c12, c12, Precision(0), Precision(0), Precision(0),
        c12, c11, c12, Precision(0), Precision(0), Precision(0),
        c12, c12, c11, Precision(0), Precision(0), Precision(0),
        Precision(0), Precision(0), Precision(0), shear, Precision(0), Precision(0),
        Precision(0), Precision(0), Precision(0), Precision(0), shear, Precision(0),
        Precision(0), Precision(0), Precision(0), Precision(0), Precision(0), shear;
    return tangent;
}

void GeneralisedIsotropicElasticity::evaluate(const AxialStrainLinearized& strain,
                                              const Precision*             state_old,
                                              Precision*                   state_new,
                                              AxialStressCauchy&           stress,
                                              Precision&                   tangent) const {
    (void) state_old;
    (void) state_new;
    tangent        = youngs;
    stress.value() = tangent * strain.value();
}

void GeneralisedIsotropicElasticity::evaluate(const AxialStrainGreenLagrange& strain,
                                              const Precision*                state_old,
                                              Precision*                      state_new,
                                              AxialStressPK2&                 stress,
                                              Precision&                      tangent) const {
    (void) state_old;
    (void) state_new;
    tangent        = youngs;
    stress.value() = tangent * strain.value();
}

void GeneralisedIsotropicElasticity::evaluate(const VolumeStrainLinearized& strain,
                                              const Precision*              state_old,
                                              Precision*                    state_new,
                                              VolumeStressCauchy&           stress,
                                              Mat6&                         tangent) const {
    (void) state_old;
    (void) state_new;
    tangent        = volume_tangent();
    stress.voigt() = tangent * strain.voigt();
}

void GeneralisedIsotropicElasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                              const Precision*                 state_old,
                                              Precision*                       state_new,
                                              VolumeStressPK2&                 stress,
                                              Mat6&                            tangent) const {
    (void) state_old;
    (void) state_new;
    tangent        = volume_tangent();
    stress.voigt() = tangent * strain.voigt();
}

void GeneralisedIsotropicElasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                              const Precision*                     state_old,
                                              Precision*                           state_new,
                                              ShellMaterialStressCauchy&            stress,
                                              Mat5&                                 tangent) const {
    (void) state_old;
    (void) state_new;
    tangent         = shell_material_tangent();
    stress.values() = tangent * strain.values();
}

void GeneralisedIsotropicElasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
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
