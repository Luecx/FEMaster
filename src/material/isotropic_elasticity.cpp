/**
 * @file isotropic_elasticity.cpp
 * @brief Implements homogeneous isotropic linear elasticity.
 *
 * The implementation constructs constant axial, plane-stress shell and
 * three-dimensional Hooke tangents. The same material operator is paired with
 * Cauchy stress for linearized kinematics and second Piola-Kirchhoff stress for
 * Green-Lagrange kinematics.
 *
 * @see IsotropicElasticity
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

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

/**
 * Constructs an isotropic Hooke material and derives its shear modulus.
 *
 * Positive Young's modulus and the open stability interval `-1 < nu < 0.5`
 * ensure finite positive shear and bulk stiffness.
 *
 * @param youngs_in Young's modulus.
 * @param poisson_in Poisson's ratio.
 */
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

/**
 * Builds the isotropic in-plane plane-stress tangent.
 *
 * The engineering ordering is `[epsilon11,epsilon22,gamma12]`; consequently the
 * shear diagonal is the engineering shear modulus rather than twice that value.
 *
 * @return Constant three-by-three plane-stress material tangent.
 */
Mat3 IsotropicElasticity::plane_stress_tangent() const {
    const Precision scalar = youngs / (Precision(1) - poisson * poisson);
    Mat3 tangent;
    tangent << scalar,           scalar * poisson, Precision(0),
               scalar * poisson, scalar,           Precision(0),
               Precision(0),     Precision(0),     shear;
    return tangent;
}

/**
 * Embeds plane-stress and transverse-shear response into shell material ordering.
 *
 * @return Constant tangent ordered `[11,22,12,13,23]`.
 */
Mat5 IsotropicElasticity::shell_material_tangent() const {
    Mat5 tangent = Mat5::Zero();
    tangent.template block<3, 3>(0, 0) = plane_stress_tangent();
    tangent(3, 3) = shear;
    tangent(4, 4) = shear;
    return tangent;
}

/**
 * Builds the full isotropic Hooke tangent in engineering-Voigt ordering.
 *
 * @return Constant six-by-six three-dimensional material tangent.
 */
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

/**
 * Evaluates linearized axial Cauchy stress from infinitesimal strain.
 *
 * @param strain Infinitesimal axial strain.
 * @param state Unused state row; isotropic Hooke elasticity is stateless.
 * @param stress Axial Cauchy stress.
 * @param tangent Constant derivative `d sigma/d epsilon = E`.
 */
void IsotropicElasticity::evaluate(const AxialStrainLinearized& strain,
                                   Precision*                   state,
                                   AxialStressCauchy&           stress,
                                   Precision&                   tangent) const {
    (void) state;
    tangent        = youngs;
    stress.value() = tangent * strain.value();
}

/**
 * Evaluates axial PK2 stress work-conjugate to Green-Lagrange strain.
 *
 * @param strain Axial Green-Lagrange strain.
 * @param state Unused state row; isotropic Hooke elasticity is stateless.
 * @param stress Axial second Piola-Kirchhoff stress.
 * @param tangent Constant material derivative `dS/dE = E`.
 */
void IsotropicElasticity::evaluate(const AxialStrainGreenLagrange& strain,
                                   Precision*                      state,
                                   AxialStressPK2&                 stress,
                                   Precision&                      tangent) const {
    (void) state;
    tangent        = youngs;
    stress.value() = tangent * strain.value();
}

/**
 * Evaluates linearized three-dimensional Cauchy stress in material coordinates.
 *
 * @param strain Infinitesimal engineering strain vector.
 * @param state Unused state row; isotropic Hooke elasticity is stateless.
 * @param stress Cauchy stress in engineering-Voigt ordering.
 * @param tangent Constant isotropic three-dimensional tangent.
 */
void IsotropicElasticity::evaluate(const VolumeStrainLinearized& strain,
                                   Precision*                    state,
                                   VolumeStressCauchy&           stress,
                                   Mat6&                         tangent) const {
    (void) state;
    tangent        = volume_tangent();
    stress.voigt() = tangent * strain.voigt();
}

/**
 * Evaluates three-dimensional PK2 stress from Green-Lagrange strain.
 *
 * @param strain Green-Lagrange engineering strain vector.
 * @param state Unused state row; isotropic Hooke elasticity is stateless.
 * @param stress Second Piola-Kirchhoff stress in material coordinates.
 * @param tangent Constant material derivative `dS/dE`.
 */
void IsotropicElasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                   Precision*                       state,
                                   VolumeStressPK2&                 stress,
                                   Mat6&                            tangent) const {
    (void) state;
    tangent        = volume_tangent();
    stress.voigt() = tangent * strain.voigt();
}

/**
 * Evaluates linearized five-component shell Cauchy stress under plane stress.
 *
 * @param strain Shell material strain ordered `[11,22,12,13,23]`.
 * @param state Unused state row; isotropic Hooke elasticity is stateless.
 * @param stress Shell Cauchy stress in the same material ordering.
 * @param tangent Constant reduced shell tangent.
 */
void IsotropicElasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                   Precision*                           state,
                                   ShellMaterialStressCauchy&            stress,
                                   Mat5&                                 tangent) const {
    (void) state;
    tangent         = shell_material_tangent();
    stress.values() = tangent * strain.values();
}

/**
 * Evaluates five-component shell PK2 stress from Green-Lagrange strain.
 *
 * @param strain Shell Green-Lagrange material strain.
 * @param state Unused state row; isotropic Hooke elasticity is stateless.
 * @param stress Shell second Piola-Kirchhoff stress.
 * @param tangent Constant reduced material derivative.
 */
void IsotropicElasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                                   Precision*                              state,
                                   ShellMaterialStressPK2&                 stress,
                                   Mat5&                                   tangent) const {
    (void) state;
    tangent         = shell_material_tangent();
    stress.values() = tangent * strain.values();
}

} // namespace fem::material
