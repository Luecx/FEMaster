/**
 * @file s3.h
 * @brief Declares the three-node linear triangular shell element.
 *
 * S3 uses the common simple-shell formulation with linear triangular geometry
 * and cubic triangular quadrature. Boundary extraction returns either the
 * positive or negative midsurface orientation by reversing the nodal winding.
 */

#pragma once

#include "shell_simple.h"
#include "../geometry/surface/surface3.h"

namespace fem::model {

struct S3 : DefaultShellElement<3,
                                Surface3,
                                math::quadrature::Domain::DOMAIN_ISO_TRI,
                                math::quadrature::Order::ORDER_CUBIC> {
    S3(ID p_elem_id, std::array<ID, 3> p_node)
        : DefaultShellElement(p_elem_id, p_node) {}

    // Recreate the concrete shell from persistent topology only. Section,
    // offsets and ModelData binding belong to the compiled Instance.
    ElementPtr copy() const override { return std::make_shared<S3>(elem_id, node_ids); }

    std::string type_name() const override { return "S3"; }

    std::shared_ptr<SurfaceInterface> surface(int surface_id) override {
        return std::make_shared<Surface3>(
            surface_id == 1
                ? std::array<ID, 3>{this->nodes()[0], this->nodes()[1], this->nodes()[2]}
                : std::array<ID, 3>{this->nodes()[0], this->nodes()[2], this->nodes()[1]}
        );
    }
};

} // namespace fem::model
