/**
 * @file abd_elasticity.h
 * @brief Declares an elasticity model backed by explicit shell ABD data.
 */

#pragma once

#include "elasticity.h"

namespace fem {
namespace material {

/**
 * @struct ABDElasticity
 * @brief Elasticity model for shell materials defined by ABD and transverse shear matrices.
 */
struct ABDElasticity : Elasticity {
    StaticMatrix<6, 6> abd = StaticMatrix<6, 6>::Zero();
    StaticMatrix<2, 2> shear = StaticMatrix<2, 2>::Zero();

    ABDElasticity(const StaticMatrix<6, 6>& abd_in,
                  const StaticMatrix<2, 2>& shear_in);

    StaticMatrix<3, 3> get_2d() override;
    StaticMatrix<6, 6> get_3d() override;
    StaticMatrix<6, 6> get_abd(Precision thickness) override;
    StaticMatrix<2, 2> get_shear(Precision thickness) override;
};

} // namespace material
} // namespace fem
