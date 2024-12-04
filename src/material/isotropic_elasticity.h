#pragma once

#include "elasticity.h"
#include "../core/core.h"

namespace fem {
namespace material {

struct IsotropicElasticity : Elasticity {
    Precision youngs;
    Precision poisson;
    Precision shear;

    IsotropicElasticity(Precision youngs, Precision poisson);

    fem::StaticMatrix<3, 3> get_2d() override;
    fem::StaticMatrix<6, 6> get_3d() override;
};

}    // namespace _material
}    // namespace fem