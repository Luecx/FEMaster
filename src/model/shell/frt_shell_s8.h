/**
 * @file frt_shell_s8.h
 * @brief Declares the eight-node finite-rotation MITC8 shell element.
 *
 * The element combines the common Total-Lagrangian finite-rotation shell
 * formulation with the classical MITC8 assumed in-layer and transverse-shear
 * fields. The quadratic serendipity geometry is evaluated directly on the
 * curved reference midsurface.
 *
 * @see FRTShell
 * @see Surface8
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#pragma once

#include "frt_shell.h"
#include "../geometry/surface/surface8.h"

namespace fem::model {

/**
 * @brief Eight-node quadratic finite-rotation MITC8 shell element.
 *
 * The first four nodes are the quadrilateral corners and the final four nodes
 * are the midside nodes on edges `0-1`, `1-2`, `2-3` and `3-0`. Membrane
 * strains and curvatures use the classical in-layer tensor interpolation,
 * while the two transverse-shear components use separate five-function
 * assumed fields.
 */
struct FRTShellS8 : FRTShell<8> {

    // Quintic quadrilateral rule corresponding to full 3 x 3 surface
    // integration.
    math::quadrature::Quadrature integration_scheme_;

    // Construct the MITC8 shell from one element identifier and eight nodal
    // identifiers in serendipity ordering.
    FRTShellS8(ID id, const std::array<ID, 8>& nodes);

    // Destroy the concrete shell through the common structural-element
    // interface.
    ~FRTShellS8() override = default;

    // Return the input/output element type identifier `MITC8FRT`.
    std::string type_name() const override;

    // Return one oriented quadratic serendipity surface representation.
    std::shared_ptr<SurfaceInterface> surface(int surface_id) override;

    // Return the full 3 x 3 quadrilateral area quadrature.
    const math::quadrature::Quadrature& integration_scheme() const override;

    // Return all numerical integration-point coordinates in `(r,s,t)` format.
    RowMatrix stress_strain_ip_rst() override;

    // Return all eight natural nodal coordinates in `(r,s,t)` format.
    RowMatrix stress_strain_nodal_rst() override;

    // Evaluate the eight quadratic serendipity shape functions.
    VecN shape_function(Precision r, Precision s) const override;

    // Evaluate the first natural derivatives of all serendipity shape functions.
    MatN2 shape_derivative(Precision r, Precision s) const override;

    // Return the natural coordinates of the four corners and four midside nodes.
    MatN2 node_coords_natural() const override;

    // Return the complete ordered MITC8 in-layer and transverse-shear sampling
    // coordinates.
    std::vector<Vec2> tying_point_coordinates() const override;

    // Replace compatible in-layer, curvature and transverse-shear components
    // and optional B rows by the classical MITC8 assumed fields.
    void apply_mitc_natural(
        const EvaluationData& data,
        const ReferencePoint& point,
        Vec8&                 strain_nat,
        Mat8x6N*              B_nat
    ) const override;

    // Apply the exact transpose of the complete MITC8 interpolation to
    // generalized natural resultant weights.
    void pull_back_mitc_resultants(
        const ReferencePoint& point,
        const Vec8&           assumed_weights,
        Vec8&                 compatible_weights,
        Span<Vec8>            tying_weights
    ) const override;
};

} // namespace fem::model
