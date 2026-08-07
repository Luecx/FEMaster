/**
 * @file generalised_isotropic_elasticity.cpp
 * @brief Implements generalized isotropic linear elasticity.
 *
 * Constant constitutive tangents combine isotropic normal coupling with the
 * independently prescribed engineering shear modulus. Linearized and
 * Green-Lagrange overloads differ only in their explicit work-conjugate stress
 * types because the underlying material law is linear in the supplied strain.
 *
 * @see GeneralisedIsotropicElasticity
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#include "generalised_isotropic_elasticity.h"

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
 * Constructs generalized isotropic elasticity from independent normal and
 * shear parameters.
 *
 * Young's modulus and Poisson's ratio satisfy the ordinary isotropic stability
 * interval, while the separately prescribed engineering shear modulus must be
 * positive.
 *
 * @param youngs_in Young's modulus controlling normal response.
 * @param poisson_in Poisson's ratio controlling normal coupling.
 * @param shear_in Independent engineering shear modulus.
 */
GeneralisedIsotropicElasticity::GeneralisedIsotropicElasticity(Precision youngs_in,
                                                               Precision poisson_in,
                                                               Precision shear_in)
    : youngs (youngs_in),
      poisson(poisson_in),
      shear  (shear_in) {
    logging::error(youngs > Precision(0),
                   "GENERALISED_ISOTROPIC: Young's modulus must be positive");
    logging::error(poisson > Precision(-1) && poisson < Precision(0.5),
                   "GENERALISED_ISOTROPIC: Poisson ratio must be in (-1, 0.5)");
    logging::error(shear > Precision(0),
                   "GENERALISED_ISOTROPIC: shear modulus must be positive");
}

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

/**
 * Builds the in-plane plane-stress tangent from isotropic normal coupling and
 * the independent engineering shear modulus.
 *
 * @return Constant tangent ordered `[11,22,12]`.
 */
Mat3 GeneralisedIsotropicElasticity::plane_stress_tangent() const {
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
Mat5 GeneralisedIsotropicElasticity::shell_material_tangent() const {
    Mat5 tangent = Mat5::Zero();
    tangent.template block<3, 3>(0, 0) = plane_stress_tangent();
    tangent(3, 3) = shear;
    tangent(4, 4) = shear;
    return tangent;
}

/**
 * Builds the full generalized isotropic engineering-Voigt tangent.
 *
 * The normal block is derived from Young's modulus and Poisson's ratio; all
 * three shear diagonals use the independently prescribed shear modulus.
 *
 * @return Constant six-by-six three-dimensional tangent.
 */
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

/**
 * Evaluates linearized axial Cauchy stress. The independent shear modulus does
 * not enter this one-dimensional response.
 *
 * @param strain Infinitesimal axial strain.
 * @param state Unused material-point state row.
 * @param stress Axial Cauchy stress.
 * @param tangent Constant derivative equal to Young's modulus.
 */
void GeneralisedIsotropicElasticity::evaluate(const AxialStrainLinearized& strain,
                                              Precision*                   state,
                                              AxialStressCauchy&           stress,
                                              Precision&                   tangent) const {
    (void) state;
    tangent        = youngs;
    stress.value() = tangent * strain.value();
}

/**
 * Evaluates axial PK2 stress from Green-Lagrange strain.
 *
 * @param strain Axial Green-Lagrange strain.
 * @param state Unused material-point state row.
 * @param stress Axial second Piola-Kirchhoff stress.
 * @param tangent Constant material derivative equal to Young's modulus.
 */
void GeneralisedIsotropicElasticity::evaluate(const AxialStrainGreenLagrange& strain,
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
 * @param state Unused material-point state row.
 * @param stress Cauchy stress in engineering-Voigt ordering.
 * @param tangent Constant generalized isotropic tangent.
 */
void GeneralisedIsotropicElasticity::evaluate(const VolumeStrainLinearized& strain,
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
 * @param state Unused material-point state row.
 * @param stress Second Piola-Kirchhoff stress in material coordinates.
 * @param tangent Constant material derivative `dS/dE`.
 */
void GeneralisedIsotropicElasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                                              Precision*                       state,
                                              VolumeStressPK2&                 stress,
                                              Mat6&                            tangent) const {
    (void) state;
    tangent        = volume_tangent();
    stress.voigt() = tangent * strain.voigt();
}

/**
 * Evaluates linearized shell Cauchy stress with independent shear response.
 *
 * @param strain Five-component shell material strain.
 * @param state Unused material-point state row.
 * @param stress Shell Cauchy stress in material ordering.
 * @param tangent Constant reduced shell tangent.
 */
void GeneralisedIsotropicElasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                                              Precision*                           state,
                                              ShellMaterialStressCauchy&            stress,
                                              Mat5&                                 tangent) const {
    (void) state;
    tangent         = shell_material_tangent();
    stress.values() = tangent * strain.values();
}

/**
 * Evaluates shell PK2 stress from five-component Green-Lagrange strain.
 *
 * @param strain Shell Green-Lagrange material strain.
 * @param state Unused material-point state row.
 * @param stress Shell second Piola-Kirchhoff stress.
 * @param tangent Constant reduced material derivative.
 */
void GeneralisedIsotropicElasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                                              Precision*                              state,
                                              ShellMaterialStressPK2&                 stress,
                                              Mat5&                                   tangent) const {
    (void) state;
    tangent         = shell_material_tangent();
    stress.values() = tangent * strain.values();
}

} // namespace fem::material
