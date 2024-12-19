#include "isotropic_elasticity.h"

fem::StaticMatrix<3, 3> fem::material::IsotropicElasticity::get_2d() {

    Precision scalar = youngs / (1 - poisson * poisson);

    return StaticMatrix<3, 3>(
        {{scalar          , scalar * poisson, 0},
         {scalar * poisson, scalar          , 0},
         {0               , 0               , (1 - poisson) / 2 * scalar}});
}

fem::StaticMatrix<6, 6> fem::material::IsotropicElasticity::get_3d() {
    Precision scalar = youngs / ((1 + poisson) * (1 - 2 * poisson));
    Precision mu     = (1 - 2 * poisson);

    return StaticMatrix<6, 6>({{1 - poisson, poisson, poisson, 0, 0, 0},
                               {poisson, 1 - poisson, poisson, 0, 0, 0},
                               {poisson, poisson, 1 - poisson, 0, 0, 0},
                               {0, 0, 0, mu / 2, 0, 0},
                               {0, 0, 0, 0, mu / 2, 0},
                               {0, 0, 0, 0, 0, mu / 2}})
           * scalar;
}

fem::StaticMatrix<2, 2>             fem::material::IsotropicElasticity::get_shear(Precision t) {
    Precision k = 5 / 6.0;
    Precision scalar = youngs * t * k / 2 / (1 + poisson);
    return StaticMatrix<2, 2>({{1, 0}, {0, 1}}) * scalar;
}
fem::StaticMatrix<3, 3>             fem::material::IsotropicElasticity::get_bend(Precision t) {
    Precision scalar = youngs * t*t*t / 12 / (1 - poisson * poisson);
    return StaticMatrix<3, 3>({
        {1, poisson, 0},
        {poisson, 1, 0},
        {0, 0, (1 - poisson) / 2}}) * scalar;
}

fem::material::IsotropicElasticity::IsotropicElasticity(Precision youngs, Precision poisson)
    : youngs(youngs)
    , poisson(poisson) {
    shear = youngs / (2 * (1 + poisson));
}
