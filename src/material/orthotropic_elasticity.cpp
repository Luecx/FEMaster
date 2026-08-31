/**
 * @file orthotropic_elasticity.cpp
 * @brief Implements homogeneous orthotropic linear elasticity.
 *
 * The implementation constructs the three-dimensional engineering compliance
 * directly from `E1`, `E2`, `E3`, `nu12`, `nu13`, `nu23`, `G12`, `G13` and
 * `G23`. Shell calls use the corresponding orthotropic plane-stress reduction.
 *
 * @see OrthotropicElasticity
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

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

OrthotropicElasticity::OrthotropicElasticity(Precision E1,
                                             Precision E2,
                                             Precision E3,
                                             Precision nu12,
                                             Precision nu13,
                                             Precision nu23,
                                             Precision G12,
                                             Precision G13,
                                             Precision G23)
    : E1(E1), E2(E2), E3(E3),
      nu12(nu12), nu13(nu13), nu23(nu23),
      G12(G12), G13(G13), G23(G23) {}

bool OrthotropicElasticity::supports_volume_linearized() const { return true; }
bool OrthotropicElasticity::supports_volume_green_lagrange() const { return true; }
bool OrthotropicElasticity::supports_shell_integration_linearized() const { return true; }
bool OrthotropicElasticity::supports_shell_integration_green_lagrange() const { return true; }

Mat3 OrthotropicElasticity::plane_stress_tangent() const {
    const Precision nu21  = nu12 * E2 / E1;
    const Precision denom = Precision(1) - nu12 * nu21;

    Mat3 tangent;
    tangent << E1 / denom,         nu12 * E2 / denom, Precision(0),
               nu21 * E1 / denom, E2 / denom,        Precision(0),
               Precision(0),       Precision(0),       G12;
    return tangent;
}

Mat5 OrthotropicElasticity::shell_material_tangent() const {
    Mat5 tangent = Mat5::Zero();
    tangent.template block<3, 3>(0, 0) = plane_stress_tangent();
    tangent(3, 3) = G13;
    tangent(4, 4) = G23;
    return tangent;
}

Mat6 OrthotropicElasticity::volume_tangent() const {
    Mat6 compliance;
    compliance <<
        Precision(1) / E1, -nu12 / E1,         -nu13 / E1,         Precision(0),       Precision(0),       Precision(0),
        -nu12 / E1,         Precision(1) / E2, -nu23 / E2,         Precision(0),       Precision(0),       Precision(0),
        -nu13 / E1,         -nu23 / E2,         Precision(1) / E3, Precision(0),       Precision(0),       Precision(0),
        Precision(0),       Precision(0),       Precision(0),       Precision(1) / G23, Precision(0),       Precision(0),
        Precision(0),       Precision(0),       Precision(0),       Precision(0),       Precision(1) / G13, Precision(0),
        Precision(0),       Precision(0),       Precision(0),       Precision(0),       Precision(0),       Precision(1) / G12;
    return compliance.inverse();
}

void OrthotropicElasticity::evaluate(const VolumeStrainLinearized& strain,
                                     const Precision*              state,
                                     VolumeStressCauchy&           stress,
                                     Mat6&                         tangent) const {
    (void) state;
    tangent        = volume_tangent();
    stress.voigt() = tangent * strain.voigt();
}

void OrthotropicElasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                     const Precision*                   state,
                                     VolumeStressPK2&                 stress,
                                     Mat6&                           tangent) const {
    (void) state;
    tangent        = volume_tangent();
    stress.voigt() = tangent * strain.voigt();
}

void OrthotropicElasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                     const Precision*                     state,
                                     ShellMaterialStressCauchy&            stress,
                                     Mat5&                               tangent) const {
    (void) state;
    tangent         = shell_material_tangent();
    stress.values() = tangent * strain.values();
}

void OrthotropicElasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                                     const Precision*                          state,
                                     ShellMaterialStressPK2&                 stress,
                                     Mat5&                                  tangent) const {
    (void) state;
    tangent         = shell_material_tangent();
    stress.values() = tangent * strain.values();
}

} // namespace fem::material
