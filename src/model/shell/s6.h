/**
 * @file s6.h
 * @brief Declares the six-node quadratic triangular shell element.
 *
 * S6 uses the common simple-shell formulation with quadratic triangular
 * geometry and cubic triangular integration. Boundary extraction preserves the
 * midsurface orientation by reversing corner and midside ordering together.
 */

#pragma once

#include "shell_simple.h"
#include "../geometry/surface/surface6.h"

namespace fem::model {

struct S6 : DefaultShellElement<6,
                                Surface6,
                                math::quadrature::Domain::DOMAIN_ISO_TRI,
                                math::quadrature::Order::ORDER_CUBIC> {
    S6(ID p_elem_id, std::array<ID, 6> p_node)
        : DefaultShellElement(p_elem_id, p_node) {}

    // Recreate only the persistent shell topology. Section assignment and all
    // compiled assembly state are attached afterwards by Model::compile().
    ElementPtr copy() const override { return std::make_shared<S6>(elem_id, node_ids); }

    std::string type_name() const override { return "S6"; }

    std::shared_ptr<SurfaceInterface> surface(int surface_id) override {
        return std::make_shared<Surface6>(
            surface_id == 1
                ? std::array<ID, 6>{this->nodes()[0], this->nodes()[1], this->nodes()[2],
                                    this->nodes()[3], this->nodes()[4], this->nodes()[5]}
                : std::array<ID, 6>{this->nodes()[0], this->nodes()[2], this->nodes()[1],
                                    this->nodes()[5], this->nodes()[4], this->nodes()[3]}
        );
    }
};

} // namespace fem::model
