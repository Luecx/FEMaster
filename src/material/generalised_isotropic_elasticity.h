/**
 * @file generalised_isotropic_elasticity.h
 * @brief Declares a generalised isotropic elasticity model with independent shear modulus.
 *
 * The model keeps the normal-stiffness part parameterised by Young's modulus
 * and Poisson's ratio while using an explicit shear modulus in 2D, 3D, and
 * shell shear/bending terms.
 *
 * @see src/material/generalised_isotropic_elasticity.cpp
 * @see src/material/elasticity.h
 */

#pragma once

#include "elasticity.h"

namespace fem {
namespace material {

/**
 * @struct GeneralisedIsotropicElasticity
 * @brief Generalised isotropic elasticity model with independent shear stiffness.
 */
struct GeneralisedIsotropicElasticity : Elasticity {
    Precision youngs;  ///< Young's modulus used for normal terms.
    Precision poisson; ///< Poisson's ratio used for normal coupling terms.
    Precision shear;   ///< Independently prescribed shear modulus.

    /**
     * @brief Constructs the generalized elasticity model.
     *
     * @param youngs_in Young's modulus.
     * @param poisson_in Poisson ratio.
     * @param shear_in Shear modulus.
     */
    GeneralisedIsotropicElasticity(Precision youngs_in, Precision poisson_in, Precision shear_in);

    /// Returns the plane-stress stiffness matrix with explicit in-plane shear modulus.
    StaticMatrix<3, 3> get_2d() override;

    /// Returns a 3D stiffness matrix using Lamé-like normal terms and explicit shear modulus.
    StaticMatrix<6, 6> get_3d() override;

    /// Returns the transverse shear stiffness matrix for shell elements.
    StaticMatrix<2, 2> get_shear(Precision thickness) override;

    /// Returns the shell bending stiffness matrix with explicit in-plane shear term.
    StaticMatrix<3, 3> get_bend(Precision thickness) override;
};

} // namespace material
} // namespace fem
