#pragma once

#include "elasticity.h"
#include "../core/core.h"

namespace fem {
namespace material {

struct OrthotropicElasticity : Elasticity {
    Precision Ex, Ey, Ez;         // Young's moduli
    Precision Gyz, Gzx, Gxy;      // Shear moduli
    Precision vyz, vzx, vxy;      // Poisson's ratios

    OrthotropicElasticity(Precision ex,
                          Precision ey,
                          Precision ez,
                          Precision gyz,
                          Precision gzx,
                          Precision gxy,
                          Precision vyz,
                          Precision vzx,
                          Precision vxy);

    fem::StaticMatrix<3, 3> get_2d() override;
    fem::StaticMatrix<6, 6> get_3d() override;

    StaticMatrix<2, 2>      get_shear(Precision t) override;
    StaticMatrix<3, 3>      get_bend(Precision t) override;
};

}    // namespace _material
}    // namespace fem