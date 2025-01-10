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
};
}
#endif //S6_H
