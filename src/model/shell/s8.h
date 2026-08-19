/**
 * @file s8.h
 * @brief Declares the eight-node quadratic serendipity shell element.
 *
 * S8 uses the common simple-shell formulation with quadratic serendipity
 * geometry and quintic quadrilateral integration. Boundary extraction reverses
 * both corner and midside ordering when the opposite midsurface is requested.
 */

#pragma once

#include "shell_simple.h"
#include "../geometry/surface/surface8.h"

namespace fem::model {

struct S8 : DefaultShellElement<8,
                                Surface8,
                                math::quadrature::Domain::DOMAIN_ISO_QUAD,
                                math::quadrature::Order::ORDER_QUINTIC> {
    S8(ID p_elem_id, std::array<ID, 8> p_node)
        : DefaultShellElement(p_elem_id, p_node) {}

    // Recreate only the persistent shell topology. Section assignment and all
    // compiled assembly state are attached afterwards by Model::compile().
    ElementPtr copy() const override { return std::make_shared<S8>(elem_id, node_ids); }

    std::string type_name() const override { return "S8"; }

    std::shared_ptr<SurfaceInterface> surface(int surface_id) override {
        return std::make_shared<Surface8>(
            surface_id == 1
                ? std::array<ID, 8>{this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3],
                                    this->nodes()[4], this->nodes()[5], this->nodes()[6], this->nodes()[7]}
                : std::array<ID, 8>{this->nodes()[0], this->nodes()[3], this->nodes()[2], this->nodes()[1],
                                    this->nodes()[7], this->nodes()[6], this->nodes()[5], this->nodes()[4]}
        );
    }
};

} // namespace fem::model
