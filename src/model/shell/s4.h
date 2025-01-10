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

};



}



#endif //S4_H
