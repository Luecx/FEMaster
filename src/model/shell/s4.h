/**
 * @file s4.h
 * @brief Declares the four-node bilinear quadrilateral shell element.
 *
 * S4 uses the common simple-shell formulation with bilinear quadrilateral
 * geometry and cubic quadrilateral integration. Boundary extraction preserves
 * the requested midsurface orientation by reversing nodal winding for the
 * negative face.
 */

#pragma once

#include "shell_simple.h"
#include "../geometry/surface/surface4.h"

namespace fem::model {

struct S4 : DefaultShellElement<4,
                                Surface4,
                                math::quadrature::Domain::DOMAIN_ISO_QUAD,
                                math::quadrature::Order::ORDER_CUBIC> {
    S4(ID p_elem_id, std::array<ID, 4> p_node)
        : DefaultShellElement(p_elem_id, p_node) {}

    // Recreate only the persistent shell topology. Section assignment and all
    // compiled assembly state are attached after cloning by Model::compile().
    ElementPtr copy() const override { return std::make_shared<S4>(elem_id, node_ids); }

    std::string type_name() const override { return "S4"; }

    std::shared_ptr<SurfaceInterface> surface(int surface_id) override {
        return std::make_shared<Surface4>(
            surface_id == 1
                ? std::array<ID, 4>{this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3]}
                : std::array<ID, 4>{this->nodes()[3], this->nodes()[2], this->nodes()[1], this->nodes()[0]}
        );
    }

    StaticMatrix<2, 12> strain_disp_shear(ShapeFunction& shape_func,
                                          ShapeDerivative& shape_der,
                                          Jacobian& jacobian) override {
        return DefaultShellElement::strain_disp_shear(shape_func, shape_der, jacobian);
    }
};

} // namespace fem::model
