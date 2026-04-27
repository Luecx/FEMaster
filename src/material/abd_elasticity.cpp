#include "abd_elasticity.h"

namespace fem {
namespace material {

ABDElasticity::ABDElasticity(const StaticMatrix<6, 6>& abd_in,
                             const StaticMatrix<2, 2>& shear_in)
    : abd(abd_in)
    , shear(shear_in) {}

StaticMatrix<3, 3> ABDElasticity::get_2d() {
    return abd.template block<3, 3>(0, 0);
}

StaticMatrix<6, 6> ABDElasticity::get_3d() {
    return StaticMatrix<6, 6>::Zero();
}

StaticMatrix<6, 6> ABDElasticity::get_abd(Precision thickness) {
    (void) thickness;
    return abd;
}

StaticMatrix<2, 2> ABDElasticity::get_shear(Precision thickness) {
    (void) thickness;
    return shear;
}

} // namespace material
} // namespace fem
