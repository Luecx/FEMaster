//
// Created by f_eggers on 10.01.2025.
//

#ifndef S8_H
#define S8_H


#include "shell_simple.h"
#include "../geometry/surface/surface8.h"
namespace fem::model {
struct S8 : DefaultShellElement<8, Surface8, quadrature::Domain::DOMAIN_ISO_QUAD, quadrature::Order::ORDER_QUINTIC> {
    S8(ID p_elem_id, std::array<ID, 8> p_node)
        : DefaultShellElement(p_elem_id, p_node) {}

    virtual std::shared_ptr<SurfaceInterface> surface(int surface_id) {
        return std::make_shared<Surface8>(
            surface_id == 1
                ? std::array<ID, 8>{this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3],
                                    this->nodes()[4], this->nodes()[5], this->nodes()[6], this->nodes()[7]}
                : std::array<ID, 8>{this->nodes()[0], this->nodes()[3], this->nodes()[2], this->nodes()[1],
                                    this->nodes()[7], this->nodes()[6], this->nodes()[5], this->nodes()[4]});
    }

};
}
#endif //S8_H
