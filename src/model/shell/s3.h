//
// Created by f_eggers on 10.01.2025.
//

#ifndef S3_H
#define S3_H

#include "shell_simple.h"
#include "../geometry/surface/surface3.h"
namespace fem::model {
struct S3 : DefaultShellElement<3, Surface3, quadrature::Domain::DOMAIN_ISO_TRI, quadrature::Order::ORDER_CUBIC> {
    S3(ID p_elem_id, std::array<ID, 3> p_node)
        : DefaultShellElement(p_elem_id, p_node) {}
};
}
#endif //S3_H
