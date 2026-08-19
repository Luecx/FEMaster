/**
 * @file line_interface.h
 * @brief Declares the non-templated interface for isoparametric line geometries.
 */

#pragma once

#include "../../../core/types_cls.h"
#include "../../../core/types_eig.h"
#include "../../../data/field.h"

namespace fem::model {

enum IsoParametricLineRange {
    MINUS_ONE_TO_ONE,
    ZERO_TO_ONE
};

struct LineInterface {
    const Index n_nodes = 0;

    explicit LineInterface(Index n)
        : n_nodes(n) {}

    virtual ~LineInterface() = default;

    virtual DynamicVector shape_function(Precision r) const = 0;
    virtual Vec3 local_to_global(Precision local, const Field& node_coords_system) const = 0;
    virtual Precision global_to_local(const Vec3& global,
                                      const Field& node_coords_system,
                                      bool clip = false) const = 0;

    virtual const ID* nodes() const = 0;
    virtual ID* nodes() = 0;
    virtual LinePtr copy() const = 0;

    const ID* begin() const { return nodes(); }
    const ID* end() const { return nodes() + n_nodes; }
    ID* begin() { return nodes(); }
    ID* end() { return nodes() + n_nodes; }
};

} // namespace fem::model
