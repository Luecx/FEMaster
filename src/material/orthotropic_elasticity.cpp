/**
 * @file orthotropic_elasticity.cpp
 * @brief Implements homogeneous orthotropic linear elasticity.
 *
 * The implementation constructs the three-dimensional engineering compliance
 * from reciprocal Poisson ratios and inverts it to obtain the material tangent.
 * Shell calls use an orthotropic plane-stress reduction and the prescribed
 * transverse shear moduli.
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
 * Constructs orthotropic elasticity from local engineering constants.
 *
 * The supplied major Poisson ratios are stored with the directional Young's and
 * shear moduli. Reciprocal minor ratios are derived when the symmetric volume
 * compliance is assembled.
 *
 * @param ex Young's modulus along material axis X.
 * @param ey Young's modulus along material axis Y.
 * @param ez Young's modulus along material axis Z.
 * @param gyz Engineering shear modulus in the YZ plane.
 * @param gzx Engineering shear modulus in the ZX plane.
 * @param gxy Engineering shear modulus in the XY plane.
 * @param vyz Major Poisson ratio for Y loading and Z contraction.
 * @param vzx Major Poisson ratio for Z loading and X contraction.
 * @param vxy Major Poisson ratio for X loading and Y contraction.
 */
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

/**
 * Builds the orthotropic in-plane plane-stress tangent.
 *
 * The reciprocal ratio `nu_yx = nu_xy E_y/E_x` enforces symmetry of the normal
 * block; the engineering XY shear term is uncoupled.
 *
 * @return Constant tangent ordered `[11,22,12]`.
 */
Mat3 OrthotropicElasticity::plane_stress_tangent() const {
    const Precision vyx   = vxy * Ey / Ex;
    const Precision denom = Precision(1) - vxy * vyx;

    Mat3 tangent;
    tangent << Ex / denom,       Ex * vyx / denom, Precision(0),
               Ey * vxy / denom, Ey / denom,       Precision(0),
               Precision(0),     Precision(0),      Gxy;
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
    tangent(3, 3) = Gzx;
    tangent(4, 4) = Gyz;
    return tangent;
}

/**
 * Constructs the full three-dimensional orthotropic tangent.
 *
 * Reciprocal Poisson ratios complete a symmetric engineering compliance matrix,
 * which is inverted to obtain the stress-strain operator in local material
 * coordinates.
 *
 * @return Constant six-by-six engineering-Voigt tangent.
 */
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
