//
// Created by f_eggers on 10.01.2025.
//

#ifndef S4_H
#define S4_H

#include "shell_simple.h"

#include "../geometry/surface/surface4.h"

namespace fem::model {
struct S4 : DefaultShellElement<4, Surface4, quadrature::Domain::DOMAIN_ISO_QUAD, quadrature::Order::ORDER_CUBIC> {
    S4(ID p_elem_id, std::array<ID, 4> p_node)
        : DefaultShellElement(p_elem_id, p_node) {}

    std::string type_name() const override { return "S4"; }

    virtual std::shared_ptr<SurfaceInterface> surface(int surface_id) override {
        return std::make_shared<Surface4>(
            surface_id == 1
                ? std::array<ID, 4> {this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3]}
                : std::array<ID, 4> {this->nodes()[3], this->nodes()[2], this->nodes()[1], this->nodes()[0]});
    }

    StaticMatrix<2, 3 * 4> strain_disp_shear(ShapeFunction& shape_func, ShapeDerivative& shape_der, Jacobian& jacobian) override {
        return DefaultShellElement::strain_disp_shear(shape_func, shape_der, jacobian);
    }
};



}



#endif //S4_H
