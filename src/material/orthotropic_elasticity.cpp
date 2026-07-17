#include "orthotropic_elasticity.h"

#include "strain/shell_material_strain_green_lagrange.h"
#include "strain/shell_material_strain_linearized.h"
#include "strain/volume_strain_green_lagrange.h"
#include "strain/volume_strain_linearized.h"
#include "stress/shell_material_stress_cauchy.h"
#include "stress/shell_material_stress_pk2.h"
#include "stress/volume_stress_cauchy.h"
#include "stress/volume_stress_pk2.h"

#include <Eigen/LU>

namespace fem::material {

OrthotropicElasticity::OrthotropicElasticity(Precision ex,
                                             Precision ey,
                                             Precision ez,
                                             Precision gyz,
                                             Precision gzx,
                                             Precision gxy,
                                             Precision vyz_in,
                                             Precision vzx_in,
                                             Precision vxy_in)
    : Ex (ex),
      Ey (ey),
      Ez (ez),
      Gyz(gyz),
      Gzx(gzx),
      Gxy(gxy),
      vyz(vyz_in),
      vzx(vzx_in),
      vxy(vxy_in) {}

bool OrthotropicElasticity::supports_volume_linearized() const {
    return true;
}

bool OrthotropicElasticity::supports_volume_green_lagrange() const {
    return true;
}

bool OrthotropicElasticity::supports_shell_integration_linearized() const {
    return true;
}

bool OrthotropicElasticity::supports_shell_integration_green_lagrange() const {
    return true;
}

Mat3 OrthotropicElasticity::plane_stress_tangent() const {
    const Precision vyx   = vxy * Ey / Ex;
    const Precision denom = Precision(1) - vxy * vyx;

    Mat3 tangent;
    tangent << Ex / denom,       Ex * vyx / denom, Precision(0),
               Ey * vxy / denom, Ey / denom,       Precision(0),
               Precision(0),     Precision(0),      Gxy;
    return tangent;
}

Mat5 OrthotropicElasticity::shell_material_tangent() const {
    Mat5 tangent = Mat5::Zero();
    tangent.template block<3, 3>(0, 0) = plane_stress_tangent();
    tangent(3, 3) = Gzx;
    tangent(4, 4) = Gyz;
    return tangent;
}

Mat6 OrthotropicElasticity::volume_tangent() const {
    const Precision vyx = vxy * Ey / Ex;
    const Precision vzy = vyz * Ez / Ey;
    const Precision vxz = vzx * Ex / Ez;

    Mat6 compliance;
    compliance <<
        Precision(1) / Ex, -vyx / Ey,          -vzx / Ez,          Precision(0),       Precision(0),       Precision(0),
        -vxy / Ex,          Precision(1) / Ey, -vzy / Ez,          Precision(0),       Precision(0),       Precision(0),
        -vxz / Ex,          -vyz / Ey,          Precision(1) / Ez, Precision(0),       Precision(0),       Precision(0),
        Precision(0),       Precision(0),       Precision(0),       Precision(1) / Gyz, Precision(0),       Precision(0),
        Precision(0),       Precision(0),       Precision(0),       Precision(0),       Precision(1) / Gzx, Precision(0),
        Precision(0),       Precision(0),       Precision(0),       Precision(0),       Precision(0),       Precision(1) / Gxy;
    return compliance.inverse();
}

void OrthotropicElasticity::evaluate(const VolumeStrainLinearized& strain,
                                     const Precision*              state_old,
                                     Precision*                    state_new,
                                     VolumeStressCauchy&           stress,
                                     Mat6&                         tangent) const {
    (void) state_old;
    (void) state_new;
    tangent        = volume_tangent();
    stress.voigt() = tangent * strain.voigt();
}

void OrthotropicElasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                     const Precision*                 state_old,
                                     Precision*                       state_new,
                                     VolumeStressPK2&                 stress,
                                     Mat6&                            tangent) const {
    (void) state_old;
    (void) state_new;
    tangent        = volume_tangent();
    stress.voigt() = tangent * strain.voigt();
}

void OrthotropicElasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                     const Precision*                     state_old,
                                     Precision*                           state_new,
                                     ShellMaterialStressCauchy&            stress,
                                     Mat5&                                 tangent) const {
    (void) state_old;
    (void) state_new;
    tangent         = shell_material_tangent();
    stress.values() = tangent * strain.values();
}

void OrthotropicElasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
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
