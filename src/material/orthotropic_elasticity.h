/**
 * @file orthotropic_elasticity.h
 * @brief Declares the orthotropic elasticity model.
 *
 * Orthotropic materials exhibit different stiffness along three orthogonal
 * axes; this model returns the corresponding stiffness matrices.
 *
 * @see src/material/orthotropic_elasticity.cpp
 * @see src/material/elasticity.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "elasticity.h"

namespace fem {
namespace material {

/**
 * @struct OrthotropicElasticity
 * @brief Elasticity model for orthotropic materials.
 */
struct OrthotropicElasticity : Elasticity {
    Precision Ex;
    Precision Ey;
    Precision Ez;
    Precision Gyz;
    Precision Gzx;
    Precision Gxy;
    Precision vyz;
    Precision vzx;
    Precision vxy;

    /**
     * @brief Constructs the orthotropic elasticity model.
     *
     * @param ex Young's modulus in the X direction.
     * @param ey Young's modulus in the Y direction.
     * @param ez Young's modulus in the Z direction.
     * @param gyz Shear modulus in the YZ plane.
     * @param gzx Shear modulus in the ZX plane.
     * @param gxy Shear modulus in the XY plane.
     * @param vyz Poisson ratio for YZ coupling.
     * @param vzx Poisson ratio for ZX coupling.
     * @param vxy Poisson ratio for XY coupling.
     */
    OrthotropicElasticity(Precision ex,
                          Precision ey,
                          Precision ez,
                          Precision gyz,
                          Precision gzx,
                          Precision gxy,
                          Precision vyz,
                          Precision vzx,
                          Precision vxy);

    /// Returns the plane-stress stiffness matrix for orthotropic material.
    StaticMatrix<3, 3> get_2d() override;

    /// Returns the full 3D stiffness matrix for orthotropic material.
    StaticMatrix<6, 6> get_3d() override;

    /// Throws because shear stiffness is not yet implemented for orthotropic shells.
    StaticMatrix<2, 2> get_shear(Precision thickness) override;

    /// Throws because bending stiffness is not yet implemented for orthotropic shells.
    StaticMatrix<3, 3> get_bend(Precision thickness) override;
};

} // namespace material
} // namespace fem
