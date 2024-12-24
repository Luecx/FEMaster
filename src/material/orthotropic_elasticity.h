/******************************************************************************
* @file orthotropic_elasticity.h
* @brief Defines the OrthotropicElasticity class for orthotropic material behavior.
*
* The OrthotropicElasticity class calculates stiffness matrices for materials with 
* different properties along three orthogonal axes. It extends the base Elasticity 
* class and includes methods for both 2D and 3D stiffness computations.
*
* @see orthotropic_elasticity.cpp
* @note Currently, shear and bending methods are not implemented.
* @date 20.12.2024
* @author Finn Eggers
******************************************************************************/

#pragma once  // Ensures this file is included only once during compilation

#include "elasticity.h"          // Base Elasticity class
#include "../core/core.h"        // Core FEM types

namespace fem::material {

/******************************************************************************
* @class OrthotropicElasticity
* @brief Class for representing orthotropic material elasticity.
*
* This class calculates stiffness matrices for orthotropic materials in 2D and 3D.
* Orthotropic materials have unique properties along three mutually orthogonal axes.
******************************************************************************/
struct OrthotropicElasticity : Elasticity {
    Precision Ex;    ///< Young's modulus in the X direction.
    Precision Ey;    ///< Young's modulus in the Y direction.
    Precision Ez;    ///< Young's modulus in the Z direction.
    Precision Gyz;   ///< Shear modulus in the YZ plane.
    Precision Gzx;   ///< Shear modulus in the ZX plane.
    Precision Gxy;   ///< Shear modulus in the XY plane.
    Precision vyz;   ///< Poisson's ratio for YZ deformation.
    Precision vzx;   ///< Poisson's ratio for ZX deformation.
    Precision vxy;   ///< Poisson's ratio for XY deformation.

    /******************************************************************************
    * @brief Constructor for OrthotropicElasticity.
    * 
    * Initializes the material properties for orthotropic behavior.
    * 
    * @param ex Young's modulus in the X direction.
    * @param ey Young's modulus in the Y direction.
    * @param ez Young's modulus in the Z direction.
    * @param gyz Shear modulus in the YZ plane.
    * @param gzx Shear modulus in the ZX plane.
    * @param gxy Shear modulus in the XY plane.
    * @param vyz Poisson's ratio for YZ deformation.
    * @param vzx Poisson's ratio for ZX deformation.
    * @param vxy Poisson's ratio for XY deformation.
    ******************************************************************************/
    OrthotropicElasticity(Precision ex,
                          Precision ey,
                          Precision ez,
                          Precision gyz,
                          Precision gzx,
                          Precision gxy,
                          Precision vyz,
                          Precision vzx,
                          Precision vxy);

    /******************************************************************************
    * @brief Computes the 2D stiffness matrix.
    * @return A 3x3 static matrix representing 2D stiffness.
    ******************************************************************************/
    fem::StaticMatrix<3, 3> get_2d() override;

    /******************************************************************************
    * @brief Computes the 3D stiffness matrix.
    * @return A 6x6 static matrix representing 3D stiffness.
    ******************************************************************************/
    fem::StaticMatrix<6, 6> get_3d() override;

    /******************************************************************************
    * @brief Computes the shear stiffness matrix (not implemented).
    * @param t Thickness of the shell element.
    * @return A 2x2 static matrix representing shear stiffness.
    * @throws std::runtime_error Always throws since not implemented.
    ******************************************************************************/
    StaticMatrix<2, 2> get_shear(Precision t) override;

    /******************************************************************************
    * @brief Computes the bending stiffness matrix (not implemented).
    * @param t Thickness of the shell element.
    * @return A 3x3 static matrix representing bending stiffness.
    * @throws std::runtime_error Always throws since not implemented.
    ******************************************************************************/
    StaticMatrix<3, 3> get_bend(Precision t) override;
};

}    // namespace fem::material
