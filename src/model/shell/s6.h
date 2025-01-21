//
// Created by f_eggers on 10.01.2025.
//

#ifndef S6_H
#define S6_H

#include "shell_simple.h"
#include "../geometry/surface/surface6.h"

namespace fem::model {
struct S6 : DefaultShellElement<6, Surface6, quadrature::Domain::DOMAIN_ISO_TRI, quadrature::Order::ORDER_CUBIC> {
    S6(ID p_elem_id, std::array<ID, 6> p_node)
        : DefaultShellElement(p_elem_id, p_node) {}

    virtual std::shared_ptr<SurfaceInterface> surface(int surface_id) {
        return std::make_shared<Surface6>(
            surface_id == 1
                ? std::array<ID, 6>{this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3],
                                    this->nodes()[4], this->nodes()[5]}
                : std::array<ID, 6>{this->nodes()[0], this->nodes()[2], this->nodes()[1], this->nodes()[5],
                                    this->nodes()[4], this->nodes()[3]});
    }

};
}
#endif //S6_H
