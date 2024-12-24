/******************************************************************************
* @file orthotropic_elasticity.cpp
* @brief Implements the OrthotropicElasticity class for orthotropic material behavior.
*
* This file contains the implementation of methods for calculating stiffness
* matrices for orthotropic materials in 2D and 3D contexts.
*
* @see orthotropic_elasticity.h
******************************************************************************/

#include "orthotropic_elasticity.h"
#include <Eigen/Eigen>  // For matrix operations

namespace fem::material {

OrthotropicElasticity::OrthotropicElasticity(Precision ex,
                                             Precision ey,
                                             Precision ez,
                                             Precision gyz,
                                             Precision gzx,
                                             Precision gxy,
                                             Precision vyz,
                                             Precision vzx,
                                             Precision vxy)
    : Ex(ex), Ey(ey), Ez(ez), Gyz(gyz), Gzx(gzx), Gxy(gxy), vyz(vyz), vzx(vzx), vxy(vxy) {}

fem::StaticMatrix<3, 3> OrthotropicElasticity::get_2d() {
    auto vyx = vxy * Ey / Ex;  // Reciprocal Poisson's ratio
    auto denom = 1 - vxy * vyx;  // Denominator term in stiffness matrix

    // Assuming plane strain condition
    return StaticMatrix<3, 3>({
        {Ex / denom, Ex * vyx / denom, 0},
        {Ey * vxy / denom, Ey / denom, 0},
        {0, 0, Gxy}});
}

fem::StaticMatrix<6, 6> OrthotropicElasticity::get_3d() {
    auto vyx = vxy * Ey / Ex;
    auto vzy = vyz * Ez / Ey;
    auto vxz = vzx * Ex / Ez;

    StaticMatrix<6, 6> S;  // Compliance matrix
    S << 1 / Ex, -vyx / Ey, -vzx / Ez, 0, 0, 0,
        -vxy / Ex, 1 / Ey, -vzy / Ez, 0, 0, 0,
        -vxz / Ex, -vyz / Ey, 1 / Ez, 0, 0, 0,
        0, 0, 0, 1 / Gyz, 0, 0,
        0, 0, 0, 0, 1 / Gzx, 0,
        0, 0, 0, 0, 0, 1 / Gxy;

    // Return the inverse of the compliance matrix to get the stiffness matrix
    return S.inverse();
}

StaticMatrix<2, 2> OrthotropicElasticity::get_shear(Precision t) {
    throw std::runtime_error("Shear stiffness is not implemented for orthotropic materials.");
}

StaticMatrix<3, 3> OrthotropicElasticity::get_bend(Precision t) {
    throw std::runtime_error("Bending stiffness is not implemented for orthotropic materials.");
}

}    // namespace fem::material
