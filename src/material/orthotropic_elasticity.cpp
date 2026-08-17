/**
 * @file orthotropic_elasticity.cpp
 * @brief Implements homogeneous orthotropic linear elasticity.
 *
 * The implementation constructs the three-dimensional engineering compliance
 * directly from `E1`, `E2`, `E3`, `nu12`, `nu13`, `nu23`, `G12`, `G13` and
 * `G23`. Symmetry supplies the reciprocal Poisson ratios implicitly through the
 * off-diagonal compliance terms. The compliance is inverted to obtain the
 * volume tangent.
 *
 * Shell calls use the corresponding orthotropic plane-stress reduction and the
 * prescribed transverse shear moduli `G13` and `G23`.
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

/**
 * Constructs orthotropic elasticity from conventional engineering constants.
 *
 * The supplied Poisson ratios are the major ratios `nu12`, `nu13` and `nu23`.
 * Reciprocal ratios are not stored independently because symmetry requires
 * `nu21/E2 = nu12/E1`, `nu31/E3 = nu13/E1` and
 * `nu32/E3 = nu23/E2`.
 *
 * @param E1 Young's modulus along material direction 1.
 * @param E2 Young's modulus along material direction 2.
 * @param E3 Young's modulus along material direction 3.
 * @param nu12 Poisson ratio for loading in direction 1 and contraction in 2.
 * @param nu13 Poisson ratio for loading in direction 1 and contraction in 3.
 * @param nu23 Poisson ratio for loading in direction 2 and contraction in 3.
 * @param G12 Engineering shear modulus in the 1-2 plane.
 * @param G13 Engineering shear modulus in the 1-3 plane.
 * @param G23 Engineering shear modulus in the 2-3 plane.
 */
OrthotropicElasticity::OrthotropicElasticity(Precision E1,
                                             Precision E2,
                                             Precision E3,
                                             Precision nu12,
                                             Precision nu13,
                                             Precision nu23,
                                             Precision G12,
                                             Precision G13,
                                             Precision G23)
    : E1  (E1),
      E2  (E2),
      E3  (E3),
      nu12(nu12),
      nu13(nu13),
      nu23(nu23),
      G12 (G12),
      G13 (G13),
      G23 (G23) {}

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

/**
 * Builds the orthotropic in-plane plane-stress tangent.
 *
 * The reciprocal ratio `nu21 = nu12 E2/E1` enforces symmetry of the normal
 * block. The engineering 1-2 shear term remains uncoupled and uses `G12`.
 *
 * @return Constant tangent ordered `[11,22,12]`.
 */
Mat3 OrthotropicElasticity::plane_stress_tangent() const {
    const Precision nu21  = nu12 * E2 / E1;
    const Precision denom = Precision(1) - nu12 * nu21;

    Mat3 tangent;
    tangent << E1 / denom,         nu12 * E2 / denom, Precision(0),
               nu21 * E1 / denom, E2 / denom,        Precision(0),
               Precision(0),       Precision(0),       G12;
    return tangent;
}

/**
 * Embeds the plane-stress block and directional transverse shear moduli into
 * the five-component shell material tangent.
 *
 * @return Tangent ordered `[11,22,12,13,23]`.
 */
Mat5 OrthotropicElasticity::shell_material_tangent() const {
    Mat5 tangent = Mat5::Zero();
    tangent.template block<3, 3>(0, 0) = plane_stress_tangent();
    tangent(3, 3) = G13;
    tangent(4, 4) = G23;
    return tangent;
}

/**
 * Constructs the full three-dimensional orthotropic tangent.
 *
 * The engineering compliance is assembled directly from the three major
 * Poisson ratios. Its off-diagonal terms are written in their symmetric form,
 * so the reciprocal minor ratios never become independent stored parameters.
 * FEMaster's volume Voigt ordering places shear components as `[23,13,12]`.
 *
 * @return Constant six-by-six engineering-Voigt tangent.
 */
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

/**
 * Evaluates linearized orthotropic Cauchy stress in material coordinates.
 *
 * @param strain Infinitesimal engineering strain vector.
 * @param state Unused material-point state row.
 * @param stress Cauchy stress in engineering-Voigt ordering.
 * @param tangent Constant orthotropic volume tangent.
 */
void OrthotropicElasticity::evaluate(const VolumeStrainLinearized& strain,
                                     Precision*                    state,
                                     VolumeStressCauchy&           stress,
                                     Mat6&                         tangent) const {
    (void) state;
    tangent        = volume_tangent();
    stress.voigt() = tangent * strain.voigt();
}

/**
 * Evaluates orthotropic PK2 stress from Green-Lagrange strain.
 *
 * @param strain Green-Lagrange engineering strain vector.
 * @param state Unused material-point state row.
 * @param stress Second Piola-Kirchhoff stress in material coordinates.
 * @param tangent Constant material derivative `dS/dE`.
 */
void OrthotropicElasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                     Precision*                       state,
                                     VolumeStressPK2&                 stress,
                                     Mat6&                            tangent) const {
    (void) state;
    tangent        = volume_tangent();
    stress.voigt() = tangent * strain.voigt();
}

/**
 * Evaluates linearized orthotropic shell Cauchy stress under plane stress.
 *
 * @param strain Five-component shell material strain.
 * @param state Unused material-point state row.
 * @param stress Shell Cauchy stress in material ordering.
 * @param tangent Constant reduced orthotropic shell tangent.
 */
void OrthotropicElasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                     Precision*                           state,
                                     ShellMaterialStressCauchy&            stress,
                                     Mat5&                                 tangent) const {
    (void) state;
    tangent         = shell_material_tangent();
    stress.values() = tangent * strain.values();
}

/**
 * Evaluates orthotropic shell PK2 stress from Green-Lagrange strain.
 *
 * @param strain Five-component Green-Lagrange material strain.
 * @param state Unused material-point state row.
 * @param stress Shell second Piola-Kirchhoff stress.
 * @param tangent Constant reduced material derivative.
 */
void OrthotropicElasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                                     Precision*                              state,
                                     ShellMaterialStressPK2&                 stress,
                                     Mat5&                                   tangent) const {
    (void) state;
    tangent         = shell_material_tangent();
    stress.values() = tangent * strain.values();
}

} // namespace fem::material
