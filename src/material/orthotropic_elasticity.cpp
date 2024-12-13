#include "orthotropic_elasticity.h"


fem::StaticMatrix<3, 3> fem::material::OrthotropicElasticity::get_2d() {
    auto vyx = vxy * Ey / Ex;  // Reciprocal Poisson's ratio

    auto denom = 1 - vxy * vyx;  // Denominator term in stiffness matrix

    // Assuming plane strain condition
    return StaticMatrix<3, 3>(
        {{
            Ex / denom, Ex * vyx / denom, 0,
            Ey * vxy / denom, Ey / denom, 0,
            0, 0, Gxy
        }}
    );
}

fem::StaticMatrix<6, 6> fem::material::OrthotropicElasticity::get_3d() {
    auto vyx = vxy * Ey / Ex;
    auto vzy = vyz * Ez / Ey;
    auto vxz = vzx * Ex / Ez;

    StaticMatrix<6,6> S;
    S << 1 / Ex, -vyx / Ey, -vzx / Ez, 0, 0, 0,
        -vxy / Ex, 1 / Ey, -vzy / Ez, 0, 0, 0,
        -vxz / Ex, -vyz / Ey, 1 / Ez, 0, 0, 0,
        0, 0, 0, 1 / Gyz, 0, 0,
        0, 0, 0, 0, 1 / Gzx, 0,
        0, 0, 0, 0, 0, 1 / Gxy;

    return S.inverse();
}

fem::StaticMatrix<2, 2> fem::material::OrthotropicElasticity::get_shear(Precision t) {
    throw std::runtime_error("Not implemented");
    return StaticMatrix<2, 2>();
}

fem::StaticMatrix<3, 3> fem::material::OrthotropicElasticity::get_bend(Precision t) {
    throw std::runtime_error("Not implemented");
    return StaticMatrix<3, 3>();
}

fem::material::OrthotropicElasticity::OrthotropicElasticity(Precision ex,
                                                            Precision ey,
                                                            Precision ez,
                                                            Precision gyz,
                                                            Precision gzx,
                                                            Precision gxy,
                                                            Precision vyz,
                                                            Precision vzx,
                                                            Precision vxy)
    : Ex(ex)
    , Ey(ey)
    , Ez(ez)
    , Gyz(gyz)
    , Gzx(gzx)
    , Gxy(gxy)
    , vyz(vyz)
    , vzx(vzx)
    , vxy(vxy) {}
