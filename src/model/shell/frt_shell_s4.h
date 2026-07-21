/**
 * @file frt_shell_s4.h
 * @brief Declares the four-node finite-rotation MITC4 shell element.
 *
 * The element combines bilinear quadrilateral geometry with the classical
 * MITC4 transverse-shear interpolation. The reference midsurface is evaluated
 * pointwise and may therefore be warped; no planar projection is introduced.
 *
 * @see FRTShell
 * @see Surface4
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#pragma once

#include "frt_shell.h"
#include "../geometry/surface/surface4.h"

namespace fem::model {

/**
 * @brief Four-node geometrically nonlinear MITC4 shell element.
 *
 * `gamma_r3` is tied at the bottom and top edge midpoints, while `gamma_s3` is
 * tied at the left and right edge midpoints. Both covariant shear components
 * are linearly interpolated and transformed into the pointwise orthonormal
 * reference basis by the common shell kernel.
 */
struct FRTShellS4 : FRTShell<4> {

    // Cubic quadrilateral area quadrature used by the shell integrations.
    math::quadrature::Quadrature integration_scheme_;

    // Construct the MITC4 shell from one element identifier and four nodal
    // identifiers in positive surface ordering.
    FRTShellS4(ID id, const std::array<ID, 4>& nodes);

    // Destroy the concrete shell through the common structural-element
    // interface.
    ~FRTShellS4() override = default;

    // Return the input/output element type identifier `MITC4FRT`.
    std::string type_name() const override;

    // Return one oriented bilinear surface representation.
    std::shared_ptr<SurfaceInterface> surface(int surface_id) override;

    // Return the cubic quadrilateral area quadrature used by this element.
    const math::quadrature::Quadrature& integration_scheme() const override;

    // Return all numerical integration-point coordinates in `(r,s,t)` format.
    RowMatrix stress_strain_ip_rst() override;

    // Return the four natural corner coordinates in `(r,s,t)` format.
    RowMatrix stress_strain_nodal_rst() override;

    // Evaluate the four bilinear quadrilateral shape functions.
    VecN shape_function(Precision r, Precision s) const override;

    // Evaluate the first natural derivatives of all bilinear shape functions.
    MatN2 shape_derivative(Precision r, Precision s) const override;

    // Return the natural coordinates of the four quadrilateral corners.
    MatN2 node_coords_natural() const override;

    // Return the four classical MITC4 edge-midpoint tying coordinates.
    std::vector<Vec2> tying_point_coordinates() const override;

    // Replace the compatible natural transverse-shear components and optional
    // B rows by the classical MITC4 assumed field.
    void apply_mitc_natural(
        const EvaluationData& data,
        const ReferencePoint& point,
        Vec8&                 strain_nat,
        Mat8x6N*              B_nat
    ) const override;

    // Apply the transpose of the MITC4 shear interpolation to natural
    // generalized resultant weights.
    void pull_back_mitc_resultants(
        const ReferencePoint& point,
        const Vec8&           assumed_weights,
        Vec8&                 compatible_weights,
        Span<Vec8>            tying_weights
    ) const override;
};

} // namespace fem::model
