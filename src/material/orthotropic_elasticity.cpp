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

    auto delta = (1 - vxy * vyx - vyz * vzy - vzx * vxz - 2 * vxy * vyz * vzx) / (Ex * Ey * Ez);

    return StaticMatrix<6, 6>(
        {{
            Ex * (  1 - vyz * vzy) * delta, Ex * (vyx + vzx * vyz) * delta, Ex * (vzx + vyx * vzy) * delta, 0  , 0  , 0,
            Ey * (vxy + vxz * vzy) * delta, Ey * (  1 - vzx * vxz) * delta, Ey * (vzy + vzx * vxy) * delta, 0  , 0  , 0,
            Ez * (vxz + vxy * vyz) * delta, Ez * (vyz + vxz * vyx) * delta, Ez * (  1 - vxy * vyx) * delta, 0  , 0  , 0,
                                         0,                              0,                              0, Gyz, 0  , 0,
                                         0,                              0,                              0, 0  , Gzx, 0,
                                         0,                              0,                              0, 0  , 0  , Gxy
        }}
    );
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
