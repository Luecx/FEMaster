/**
 * @file frt_shell_s6.h
 * @brief Declares the six-node finite-rotation MITC6-b shell element.
 *
 * The element uses a quadratic isoparametric triangular midsurface and the
 * MITC6-b assumed-strain interpolation. Both in-plane and transverse-shear
 * natural strain components are mixed interpolated to suppress membrane and
 * shear locking. The same interpolation acts on membrane strain, curvature, B
 * rows and transposed geometric-resultant weights.
 *
 * @see FRTShell
 * @see Surface6
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#pragma once

#include "frt_shell.h"
#include "../geometry/surface/surface6.h"

namespace fem::model {

/**
 * @brief Six-node geometrically nonlinear MITC6-b shell element.
 *
 * The quadratic triangle follows the isotropic MITC6 in-plane interpolation
 * and the linear MITC6-b transverse-shear field. Two Gauss tying positions are
 * used on every edge together with additional in-plane sampling points.
 */
struct FRTShellS6 : FRTShell<6> {

    // Cubic triangular area quadrature used by all shell integrations.
    math::quadrature::Quadrature integration_scheme_;

    // Construct the MITC6-b shell from one element identifier and six nodal
    // identifiers in quadratic triangular ordering.
    FRTShellS6(ID id, const std::array<ID, 6>& nodes);

    // Destroy the concrete shell through the common structural-element
    // interface.
    ~FRTShellS6() override = default;

    // Return the input/output element type identifier `MITC6FRT`.
    std::string type_name() const override;

    // Return one oriented quadratic triangular surface representation.
    std::shared_ptr<SurfaceInterface> surface(int surface_id) override;

    // Return the cubic triangular area quadrature used by this element.
    const math::quadrature::Quadrature& integration_scheme() const override;

    // Return all numerical integration-point coordinates in `(r,s,t)` format.
    RowMatrix stress_strain_ip_rst() override;

    // Return all six natural nodal coordinates in `(r,s,t)` format.
    RowMatrix stress_strain_nodal_rst() override;

    // Evaluate the six quadratic triangular shape functions.
    VecN shape_function(Precision r, Precision s) const override;

    // Evaluate the first natural derivatives of all quadratic shape functions.
    MatN2 shape_derivative(Precision r, Precision s) const override;

    // Return the natural coordinates of the three corners and three midside
    // nodes.
    MatN2 node_coords_natural() const override;

    // Return the ordered MITC6 in-plane and MITC6-b transverse-shear sampling
    // coordinates.
    std::vector<Vec2> tying_point_coordinates() const override;

    // Replace compatible in-plane, curvature and transverse-shear components
    // and optional B rows by the MITC6-b assumed fields.
    void apply_mitc_natural(
        const EvaluationData& data,
        const ReferencePoint& point,
        Vec8&                 strain_nat,
        Mat8x6N*              B_nat
    ) const override;

    // Apply the exact transpose of the complete MITC6-b interpolation to
    // generalized natural resultant weights.
    void pull_back_mitc_resultants(
        const ReferencePoint& point,
        const Vec8&           assumed_weights,
        Vec8&                 compatible_weights,
        Span<Vec8>            tying_weights
    ) const override;
};

} // namespace fem::model
