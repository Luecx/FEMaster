/**
 * @file frt_shell_s3.h
 * @brief Declares the three-node finite-rotation MITC3 shell element.
 *
 * The element combines the common Total-Lagrangian finite-rotation shell kernel
 * with the isotropic MITC3 transverse-shear interpolation. The geometry is a
 * linear triangular midsurface, and the three tying points lie at the edge
 * midpoints of the natural triangle.
 *
 * The implementation uses the original MITC3 assumed shear field and does not
 * introduce the internal cubic bubble rotations of MITC3+.
 *
 * @see FRTShell
 * @see Surface3
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#pragma once

#include "frt_shell.h"
#include "../geometry/surface/surface3.h"

namespace fem::model {

/**
 * @brief Three-node geometrically nonlinear MITC3 shell element.
 *
 * The element uses linear triangular geometry and an isotropic assumed
 * transverse-shear field that is constant along all three element edges. The
 * normal edge components are tied at the edge midpoints and interpolated in a
 * rotated Raviart-Thomas field over the natural triangle.
 */
struct FRTShellS3 : FRTShell<3> {

    // Cubic triangular area quadrature used for stiffness, mass and resultants.
    math::quadrature::Quadrature integration_scheme_;

    // Construct the MITC3 shell from one element identifier and three nodal
    // identifiers in counter-clockwise positive-surface ordering.
    FRTShellS3(ID id, const std::array<ID, 3>& nodes);

    // Destroy the concrete shell through the common structural-element
    // interface.
    ~FRTShellS3() override = default;

    // Return the input/output element type identifier `MITC3FRT`.
    std::string type_name() const override;

    // Return one oriented triangular surface representation for contact,
    // coupling and distributed surface operations.
    std::shared_ptr<SurfaceInterface> surface(int surface_id) override;

    // Return the cubic triangular area quadrature used by this element.
    const math::quadrature::Quadrature& integration_scheme() const override;

    // Return all numerical integration-point coordinates in `(r,s,t)` format.
    RowMatrix stress_strain_ip_rst() override;

    // Return the three natural nodal coordinates in `(r,s,t)` format.
    RowMatrix stress_strain_nodal_rst() override;

    // Evaluate the three linear barycentric shape functions.
    VecN shape_function(Precision r, Precision s) const override;

    // Evaluate the constant natural derivatives of the linear shape functions.
    MatN2 shape_derivative(Precision r, Precision s) const override;

    // Return the natural coordinates of the three triangle vertices.
    MatN2 node_coords_natural() const override;

    // Return the three edge-midpoint tying coordinates of the MITC3 shear field.
    std::vector<Vec2> tying_point_coordinates() const override;

    // Replace the two compatible natural transverse-shear components and their
    // optional B rows by the isotropic MITC3 assumed field.
    void apply_mitc_natural(
        const EvaluationData& data,
        const ReferencePoint& point,
        Vec8&                 strain_nat,
        Mat8x6N*              B_nat
    ) const override;

    // Apply the transpose of the MITC3 shear interpolation to generalized
    // natural resultant weights without allocating temporary storage.
    void pull_back_mitc_resultants(
        const ReferencePoint& point,
        const Vec8&           assumed_weights,
        Vec8&                 compatible_weights,
        Span<Vec8>            tying_weights
    ) const override;
};

} // namespace fem::model
