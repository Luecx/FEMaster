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
#include "strain/axial_strain.h"
#include "strain/shell_material_strain.h"
#include "strain/volume_strain.h"
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
 * Evaluates axial PK2 stress work-conjugate to Green-Lagrange strain.
 *
 * The same linear material law is interpreted in the reference configuration,
 *
 *     S = E E_GL
 *     dS/dE_GL = E.
 *
 * @param strain Axial Green-Lagrange strain.
 * @param old_state Unused input state row; isotropic Hooke elasticity is stateless.
 * @param new_state Unused output state row.
 * @param stress Axial second Piola-Kirchhoff stress.
 * @param tangent Optional material derivative `dS/dE`.
 */
void IsotropicElasticity::evaluate(const AxialStrain&              strain,
                                   const Precision*                old_state,
                                   Precision*                      new_state,
                                   AxialStressPK2&                 stress,
                                   Precision*                      tangent) const {
    (void) old_state;
    (void) new_state;

    // Evaluate PK2 stress directly from the supplied Green-Lagrange strain.
    stress.value() = youngs * strain.value();

    if (tangent != nullptr) {
        *tangent = youngs;
    }
}

/**
 * Evaluates three-dimensional PK2 stress from Green-Lagrange strain.
 *
 * @param strain Green-Lagrange engineering strain vector.
 * @param old_state Unused input state row; isotropic Hooke elasticity is stateless.
 * @param new_state Unused output state row.
 * @param stress Second Piola-Kirchhoff stress in material coordinates.
 * @param tangent Optional material derivative `dS/dE`.
 */
void IsotropicElasticity::evaluate(const VolumeStrain&              strain,
                                   const Precision*                 old_state,
                                   Precision*                       new_state,
                                   VolumeStressPK2&                 stress,
                                   Mat6*                            tangent) const {
    (void) old_state;
    (void) new_state;

    // The Saint-Venant-Kirchhoff response uses the same constant Hooke operator.
    const Mat6 material_tangent = volume_tangent();
    stress.voigt() = material_tangent * strain.voigt();

    if (tangent != nullptr) {
        *tangent = material_tangent;
    }
}

/**
 * Evaluates five-component shell PK2 stress from Green-Lagrange strain.
 *
 * @param strain Shell Green-Lagrange material strain.
 * @param old_state Unused input state row; isotropic Hooke elasticity is stateless.
 * @param new_state Unused output state row.
 * @param stress Shell second Piola-Kirchhoff stress.
 * @param tangent Optional reduced material derivative.
 */
void IsotropicElasticity::evaluate(const ShellMaterialStrain&              strain,
                                   const Precision*                        old_state,
                                   Precision*                              new_state,
                                   ShellMaterialStressPK2&                 stress,
                                   Mat5*                                   tangent) const {
    (void) old_state;
    (void) new_state;

    const Mat5 material_tangent = shell_material_tangent();
    stress.values() = material_tangent * strain.values();

    if (tangent != nullptr) {
        *tangent = material_tangent;
    }
}

} // namespace fem::material
