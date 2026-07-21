/**
 * @file frt_shell_s4.h
 * @brief Declares the four-node finite-rotation MITC4 shell element.
 *
 * The element combines bilinear quadrilateral geometry with the classical
 * MITC4 transverse-shear interpolation. The reference midsurface is evaluated
 * pointwise and may therefore be warped; no planar element projection is used.
 *
 * @see FRTShell
 * @see Surface4
 *
 * @author Finn Eggers
 * @date 20.07.2026
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
 * are linearly interpolated and transformed into the local pointwise reference
 * basis by the common shell kernel.
 */
struct FRTShellS4 : FRTShell<4> {
    // Element geometry and quadrilateral area quadrature
    Surface4                     geometry;
    math::quadrature::Quadrature integration_scheme_;

    // Construction
    FRTShellS4(ID id, const std::array<ID, 4>& nodes);
    ~FRTShellS4() override = default;

    // Element identification, surface extraction and numerical integration
    std::string type_name() const override;
    std::shared_ptr<SurfaceInterface> surface(int surface_id) override;
    const math::quadrature::Quadrature& integration_scheme() const override;

    // Integration-point and nodal output coordinates
    RowMatrix stress_strain_ip_rst() override;
    RowMatrix stress_strain_nodal_rst() override;

    // Bilinear interpolation and MITC4 tying
    VecN  shape_function      (Precision r, Precision s) const override;
    MatN2 shape_derivative    (Precision r, Precision s) const override;
    MatN2 node_coords_natural () const override;
    std::vector<Vec2> tying_point_coordinates() const override;

    void apply_mitc_natural(const EvaluationData& data,
                            const ReferencePoint& point,
                            Vec8&                 strain_nat,
                            Mat8x6N*              B_nat,
                            Vec6NMat*             G_nat) const override;
};

} // namespace fem::model
