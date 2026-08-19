/**
 * @file c3d8.h
 * @brief Declares the eight-node trilinear hexahedral solid element.
 *
 * `C3D8` supplies the natural-coordinate interpolation, derivatives, nodal
 * coordinates, boundary-face extraction and volume integration scheme required
 * by the generic `SolidElement<8>` implementation. Constitutive response,
 * nonlinear kinematics and element-level structural assembly remain implemented
 * by the solid-element base classes.
 *
 * Boundary surfaces follow the element's six face ids and are returned with
 * outward node winding so their geometric normals are suitable for pressure,
 * tie and surface-to-surface contact formulations.
 *
 * @see C3D8
 * @see SolidElement
 * @see Surface4
 *
 * @author Finn Eggers
 * @date 10.08.2026
 */

#pragma once

#include "element_solid.h"

namespace fem {
namespace model {

/**
 * @brief Eight-node trilinear isoparametric hexahedral solid element.
 *
 * The element uses the natural cube `(r, s, t) in [-1, 1]^3` with one node at
 * each corner. `SolidElement<8>` provides the common solid mechanics while this
 * specialization defines the C3D8 interpolation and topology. Six quadrilateral
 * boundary faces are available through `surface()` and use outward orientation.
 */
struct C3D8 : public SolidElement<8> {
    // Construction and element identification
    C3D8(ID pElemId, const std::array<ID, 8>& pNodeIds);

    // Recreate the concrete element from persistent topology only. Dense ids,
    // sections, offsets and ModelData binding are instance-specific compile data.
    ElementPtr copy() const override { return std::make_shared<C3D8>(elem_id, node_ids); }

    std::string type_name() const override { return "C3D8"; }

    // Natural-coordinate interpolation and reference node positions
    StaticMatrix<8, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<8, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<8, 3> node_coords_local() override;

    // Boundary topology and element integration rule
    SurfacePtr surface(ID surface_id) override;
    const math::quadrature::Quadrature& integration_scheme() const override;
};

} // namespace model
} // namespace fem
