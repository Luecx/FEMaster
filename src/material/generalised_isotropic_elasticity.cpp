/**
 * @file generalised_isotropic_elasticity.cpp
 * @brief Implements the generalised isotropic elasticity model.
 *
 * This constitutive law keeps the normal terms driven by Young's modulus and
 * Poisson's ratio while prescribing shear stiffness independently.
 *
 * @see src/material/generalised_isotropic_elasticity.h
 */

#include "generalised_isotropic_elasticity.h"

namespace fem {
namespace material {

GeneralisedIsotropicElasticity::GeneralisedIsotropicElasticity(Precision youngs_in,
                                                               Precision poisson_in,
                                                               Precision shear_in)
    : youngs(youngs_in)
    , poisson(poisson_in)
    , shear(shear_in) {}

StaticMatrix<3, 3> GeneralisedIsotropicElasticity::get_2d() {
    const Precision scalar = youngs / (1 - poisson * poisson);
    return StaticMatrix<3, 3>({
        {scalar, scalar * poisson, 0},
        {scalar * poisson, scalar, 0},
        {0, 0, shear}});
}

StaticMatrix<6, 6> GeneralisedIsotropicElasticity::get_3d() {
    const Precision scalar = youngs / ((1 + poisson) * (1 - 2 * poisson));
    const Precision c11 = (1 - poisson) * scalar;
    const Precision c12 = poisson * scalar;
    return StaticMatrix<6, 6>({
        {c11, c12, c12, 0, 0, 0},
        {c12, c11, c12, 0, 0, 0},
        {c12, c12, c11, 0, 0, 0},
        {0, 0, 0, shear, 0, 0},
        {0, 0, 0, 0, shear, 0},
        {0, 0, 0, 0, 0, shear}});
}

} // namespace material
} // namespace fem
