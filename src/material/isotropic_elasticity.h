/******************************************************************************
* @file isotropic_elasticity.h
* @brief Defines the IsotropicElasticity class for isotropic material behavior.
*
* The IsotropicElasticity class computes the material stiffness matrices for 
* isotropic materials under 2D, 3D, shear, and bending conditions. It derives 
* from the Elasticity base class and uses material properties such as Young's 
* modulus and Poisson's ratio to calculate the stiffness matrices.
*
* @see isotropic_elasticity.cpp
* @author Finn Eggers
* @date 20.12.2024
******************************************************************************/

#pragma once  // Ensures this file is included only once during compilation

#include "elasticity.h"          // Base Elasticity class
#include "../core/core.h"        // Core FEM types

namespace fem::material {

/******************************************************************************
* @class IsotropicElasticity
* @brief Class representing isotropic material elasticity.
*
* This class provides methods to compute stiffness matrices for 2D and 3D 
* isotropic materials, as well as for shear and bending in plate and shell 
* formulations. The material behavior is defined by Young's modulus and 
* Poisson's ratio.
******************************************************************************/
struct IsotropicElasticity : Elasticity {
    Precision youngs;    ///< Young's modulus (E).
    Precision poisson;   ///< Poisson's ratio (ν).
    Precision shear;     ///< Shear modulus (G).

    /******************************************************************************
    * @brief Constructor for the IsotropicElasticity class.
    *
    * Initializes the material properties based on Young's modulus and Poisson's ratio.
    * The shear modulus is computed automatically.
    *
    * @param youngs Young's modulus (E).
    * @param poisson Poisson's ratio (ν).
    ******************************************************************************/
    IsotropicElasticity(Precision youngs, Precision poisson);

    /******************************************************************************
    * @brief Computes the 2D stiffness matrix for isotropic material.
    * @return A 3x3 static matrix representing 2D stiffness.
    ******************************************************************************/
    fem::StaticMatrix<3, 3> get_2d() override;

    /******************************************************************************
    * @brief Computes the 3D stiffness matrix for isotropic material.
    * @return A 6x6 static matrix representing 3D stiffness.
    ******************************************************************************/
    fem::StaticMatrix<6, 6> get_3d() override;

    /******************************************************************************
    * @brief Computes the shear stiffness matrix for isotropic material.
    * @param t Thickness of the material.
    * @return A 2x2 static matrix representing shear stiffness.
    ******************************************************************************/
    StaticMatrix<2, 2> get_shear(Precision t) override;

    /******************************************************************************
    * @brief Computes the bending stiffness matrix for isotropic material.
    * @param t Thickness of the material.
    * @return A 3x3 static matrix representing bending stiffness.
    ******************************************************************************/
    StaticMatrix<3, 3> get_bend(Precision t) override;
};

}    // namespace fem::material
