/**
 * @file line_interface.h
 * @brief Declares the non-templated interface for isoparametric line geometries.
 *
 * The interface exposes the operations required by generic model code without
 * depending on the interpolation order or natural-coordinate range of a line
 * element. The templated interpolation and geometry implementation is defined
 * separately in `line.h` and `line.inl`.
 *
 * @see Line
 * @see LineInterface
 */

#pragma once

#include "../../../core/types_eig.h"
#include "../../../data/field.h"

namespace fem {
namespace model {

/**
 * @enum IsoParametricLineRange
 * @brief Specifies the natural-coordinate range for isoparametric line elements.
 *
 * The selected range determines the natural integration domain and the
 * coordinate transformation used by the concrete shape functions.
 */
enum IsoParametricLineRange {
    MINUS_ONE_TO_ONE,
    ZERO_TO_ONE
};

/**
 * @brief Non-templated interface shared by all finite-element line elements.
 *
 * `LineInterface` provides generic access to shape functions, coordinate
 * transformations and contiguous node connectivity. The node count is stored
 * at runtime so model code can traverse line elements without knowing their
 * interpolation order.
 */
struct LineInterface {
    // Number of nodes used by the concrete line element. Generic model code
    // uses this runtime value to traverse connectivity without knowing the
    // interpolation order of the derived element.
    const Index n_nodes = 0;

    explicit LineInterface(Index n)
        : n_nodes(n) {}

    virtual ~LineInterface() = default;

    // Evaluate the element interpolation and map between natural and global
    // coordinates. Implementations use the element's natural-coordinate range
    // and nodal field to construct the physical line position.
    virtual DynamicVector shape_function(Precision r) const = 0;
    virtual Vec3      local_to_global(Precision local, const Field& node_coords_system) const = 0;
    virtual Precision global_to_local(
        const Vec3&  global,
        const Field& node_coords_system,
        bool         clip = false
    ) const = 0;

    // Access contiguous global node connectivity for const and mutable clients.
    // The returned storage follows the local interpolation order and is also
    // exposed through the begin()/end() range helpers below.
    virtual const ID* nodes() const = 0;
    virtual ID*       nodes()       = 0;

    // Provide iterator-style access to exactly n_nodes identifiers. These
    // helpers allow generic algorithms to inspect line connectivity without
    // depending on the concrete line type.
    const ID* begin() const { return nodes(); }
    const ID* end()   const { return nodes() + n_nodes; }
    ID* begin() { return nodes(); }
    ID* end()   { return nodes() + n_nodes; }
};

} // namespace model
} // namespace fem
