/**
 * @file isotropic_elasticity.cpp
 * @brief Implements the isotropic elasticity model.
 *
 * Provides stiffness matrices for isotropic materials across various shell and
 * solid formulations.
 *
 * @see src/material/isotropic_elasticity.h
 * @see src/material/elasticity.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "isotropic_elasticity.h"

namespace fem {
namespace material {

/**
 * @copydoc IsotropicElasticity::IsotropicElasticity
 */
IsotropicElasticity::IsotropicElasticity(Precision youngs_in, Precision poisson_in)
    : youngs(youngs_in)
    , poisson(poisson_in)
    , shear(youngs_in / (2 * (1 + poisson_in))) {}

/**
 * @copydoc IsotropicElasticity::get_2d
 */
StaticMatrix<3, 3> IsotropicElasticity::get_2d() {
    Precision scalar = youngs / (1 - poisson * poisson);
    return StaticMatrix<3, 3>({
        {scalar, scalar * poisson, 0},
        {scalar * poisson, scalar, 0},
        {0, 0, (1 - poisson) / 2 * scalar}});
}

/**
 * @copydoc IsotropicElasticity::get_3d
 */
StaticMatrix<6, 6> IsotropicElasticity::get_3d() {
    Precision scalar = youngs / ((1 + poisson) * (1 - 2 * poisson));
    Precision mu = (1 - 2 * poisson);
    return StaticMatrix<6, 6>({
               {1 - poisson, poisson, poisson, 0, 0, 0},
               {poisson, 1 - poisson, poisson, 0, 0, 0},
               {poisson, poisson, 1 - poisson, 0, 0, 0},
               {0, 0, 0, mu / 2, 0, 0},
               {0, 0, 0, 0, mu / 2, 0},
               {0, 0, 0, 0, 0, mu / 2}}) * scalar;
}

} // namespace material
} // namespace fem
