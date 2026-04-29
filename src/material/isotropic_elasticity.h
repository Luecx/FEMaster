/**
 * @file isotropic_elasticity.h
 * @brief Declares the isotropic elasticity model.
 *
 * Computes stiffness matrices for isotropic materials in 2D/3D.
 *
 * @see src/material/isotropic_elasticity.cpp
 * @see src/material/elasticity.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "elasticity.h"

namespace fem {
namespace material {

/**
 * @struct IsotropicElasticity
 * @brief Elasticity model for isotropic materials.
 */
struct IsotropicElasticity : Elasticity {
    Precision youngs;  ///< Young's modulus (E).
    Precision poisson; ///< Poisson's ratio (nu).
    Precision shear;   ///< Shear modulus (G).

    /**
     * @brief Constructs the isotropic elasticity model.
     *
     * @param youngs_in Young's modulus.
     * @param poisson_in Poisson's ratio.
     */
    IsotropicElasticity(Precision youngs_in, Precision poisson_in);

    /// Returns the plane-stress stiffness matrix.
    StaticMatrix<3, 3> get_2d() override;

    /// Returns the full 3D stiffness matrix.
    StaticMatrix<6, 6> get_3d() override;
};
} // namespace material
} // namespace fem
