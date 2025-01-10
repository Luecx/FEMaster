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


};
}
#endif //S8_H
