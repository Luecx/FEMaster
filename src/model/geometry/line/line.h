/**
 * @file line.h
 * @brief Declares the common templated interface for isoparametric line elements.
 *
 * `Line<N, CR>` stores fixed-size connectivity and declares the operations
 * shared by two-, three- and four-node isoparametric lines. Concrete line types
 * provide the fixed-size shape functions and their natural-coordinate
 * derivatives. The common definitions for interpolation, length integration,
 * closest-point projection and coordinate mapping are implemented in
 * `line.inl`.
 *
 * @tparam N Number of line nodes.
 * @tparam CR Natural-coordinate range, either `MINUS_ONE_TO_ONE` or
 *            `ZERO_TO_ONE`.
 *
 * @see LineInterface
 * @see line.inl
 */

#pragma once

#include "line_interface.h"

#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "../../../math/quadrature.h"

namespace fem {
namespace model {

/**
 * @brief Common templated interface for isoparametric line elements.
 *
 * The template parameter `N` selects the interpolation order and `CR` selects
 * the natural-coordinate interval. Derived classes must provide the fixed-size
 * shape functions and their first and second derivatives. This class exposes
 * the element connectivity and declares the generic geometry operations that
 * are defined in `line.inl`.
 *
 * @tparam N Number of nodes in the line element. Supported values are 2, 3
 *           and 4.
 * @tparam CR Natural-coordinate range of the element.
 */
template<Index N, IsoParametricLineRange CR = MINUS_ONE_TO_ONE>
struct Line : LineInterface {
    static_assert(N == 2 || N == 3 || N == 4,
        "Number of nodes (N) must be 2, 3, or 4.");

    // Compile-time number of nodes in the line element. This value determines
    // the interpolation storage and must agree with the concrete element's
    // topology and global connectivity.
    constexpr static Index num_nodes = N;

    // Global identifiers of the line nodes in element-local interpolation
    // order. Derived shape functions and all geometry operations use this same
    // ordering when gathering nodal coordinates.
    std::array<ID, N> node_ids;

    // Construction from the global identifiers of the element nodes. The
    // constructor initializes the non-templated interface with the same node
    // count exposed by num_nodes.
    explicit Line(const std::array<ID, N>& node_ids);

    // Contiguous connectivity access for generic model traversal. The pointers
    // refer directly to the fixed-size node identifier storage above.
    const ID* nodes() const override;
    ID*       nodes()       override;

    // Element-specific interpolation data. Concrete line types provide the
    // fixed-size shape functions and natural-coordinate derivatives used to
    // construct physical positions, tangents, curvatures and quadrature data.
    virtual StaticMatrix<N, 1> _shape_function(Precision r) const = 0;
    virtual StaticMatrix<N, 1> shape_derivative(Precision r) const = 0;
    virtual StaticMatrix<N, 1> shape_second_derivative(Precision r) const = 0;

    // Generic access to shape-function values and nodal coordinates. These
    // operations adapt concrete fixed-size interpolation data to the common
    // dynamic interface and gather coordinates from the global nodal field.
    DynamicVector shape_function(Precision r) const override;
    virtual StaticMatrix<N, 3> node_coords_global(const Field& node_coords_system) const;

    // Physical line geometry and coordinate transformations. The definitions
    // integrate the physical tangent over the selected natural range and solve
    // the closest-point problem by Newton iteration in that same coordinate
    // system.
    Precision length(const Field& node_coords_system) const;
    Precision global_to_local(
        const Vec3&  p,
        const Field& node_coords_system,
        bool         clip = true
    ) const override;
    Vec3 local_to_global(Precision r, const Field& node_coords_system) const override;

    // Natural-coordinate limits selected by CR. These bounds are used by
    // projection clipping and by callers that need to sample the valid domain.
    constexpr Precision min_r() const;
    constexpr Precision max_r() const;

    // Return the element-local node identifiers in their interpolation order.
    // The returned array is the authoritative connectivity representation used
    // by the line geometry and boundary-projection implementations.
    const std::array<ID, N>& get_node_ids() const;
};

} // namespace model
} // namespace fem

#include "line.inl"
