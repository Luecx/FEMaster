/******************************************************************************
 * @file orthotropic_elasticity.cpp
 * @brief Implements the orthotropic elasticity model.
 *
 * Provides stiffness matrices for orthotropic materials and notes unimplemented
 * shell-specific behaviour.
 *
 * @see src/material/orthotropic_elasticity.h
 * @see src/material/elasticity.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#include "orthotropic_elasticity.h"

#include <Eigen/Eigen>
#include <stdexcept>

namespace fem {
namespace material {

/******************************************************************************
 * @copydoc OrthotropicElasticity::OrthotropicElasticity
 ******************************************************************************/
OrthotropicElasticity::OrthotropicElasticity(Precision ex,
                                             Precision ey,
                                             Precision ez,
                                             Precision gyz,
                                             Precision gzx,
                                             Precision gxy,
                                             Precision vyz_in,
                                             Precision vzx_in,
                                             Precision vxy_in)
    : Ex(ex)
    , Ey(ey)
    , Ez(ez)
    , Gyz(gyz)
    , Gzx(gzx)
    , Gxy(gxy)
    , vyz(vyz_in)
    , vzx(vzx_in)
    , vxy(vxy_in) {}

/******************************************************************************
 * @copydoc OrthotropicElasticity::get_2d
 ******************************************************************************/
StaticMatrix<3, 3> OrthotropicElasticity::get_2d() {
    const Precision vyx = vxy * Ey / Ex;
    const Precision denom = 1 - vxy * vyx;
    return StaticMatrix<3, 3>({
        {Ex / denom, Ex * vyx / denom, 0},
        {Ey * vxy / denom, Ey / denom, 0},
        {0, 0, Gxy}});
}

/******************************************************************************
 * @copydoc OrthotropicElasticity::get_3d
 ******************************************************************************/
StaticMatrix<6, 6> OrthotropicElasticity::get_3d() {
    const Precision vyx = vxy * Ey / Ex;
    const Precision vzy = vyz * Ez / Ey;
    const Precision vxz = vzx * Ex / Ez;

    StaticMatrix<6, 6> compliance;
    compliance <<
        1 / Ex, -vyx / Ey, -vzx / Ez, 0, 0, 0,
        -vxy / Ex, 1 / Ey, -vzy / Ez, 0, 0, 0,
        -vxz / Ex, -vyz / Ey, 1 / Ez, 0, 0, 0,
        0, 0, 0, 1 / Gyz, 0, 0,
        0, 0, 0, 0, 1 / Gzx, 0,
        0, 0, 0, 0, 0, 1 / Gxy;

    return compliance.inverse();
}

/******************************************************************************
 * @copydoc OrthotropicElasticity::get_shear
 ******************************************************************************/
StaticMatrix<2, 2> OrthotropicElasticity::get_shear(Precision) {
    throw std::runtime_error("Shear stiffness is not implemented for orthotropic materials.");
}

/******************************************************************************
 * @copydoc OrthotropicElasticity::get_bend
 ******************************************************************************/
StaticMatrix<3, 3> OrthotropicElasticity::get_bend(Precision) {
    throw std::runtime_error("Bending stiffness is not implemented for orthotropic materials.");
}

} // namespace material
} // namespace fem
