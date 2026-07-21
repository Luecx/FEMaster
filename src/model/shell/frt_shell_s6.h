/**
 * @file frt_shell_s6.h
 * @brief Declares the six-node finite-rotation MITC6-b shell element.
 *
 * The element uses a quadratic isoparametric triangular midsurface and the
 * MITC6-b assumed-strain interpolation. Both in-plane and transverse-shear
 * natural strain components are mixed interpolated to suppress membrane and
 * shear locking. The same interpolation is applied independently to membrane
 * strain, curvature, B and G.
 *
 * @see FRTShell
 * @see Surface6
 *
 * @author Finn Eggers
 * @date 20.07.2026
 */

#pragma once

#include "frt_shell.h"
#include "../geometry/surface/surface6.h"

namespace fem::model {

/**
 * @brief Six-node geometrically nonlinear MITC6-b shell element.
 *
 * The quadratic triangle follows the isotropic MITC6 in-plane interpolation
 * and the linear MITC6-b transverse-shear field. The element uses two Gauss
 * tying positions on each edge and additional interior points for the in-plane
 * strain interpolation.
 */
struct FRTShellS6 : FRTShell<6> {
    // Element geometry and triangular area quadrature
    Surface6                     geometry;
    math::quadrature::Quadrature integration_scheme_;

    // Construction
    FRTShellS6(ID id, const std::array<ID, 6>& nodes);
    ~FRTShellS6() override = default;

    // Element identification, surface extraction and numerical integration
    std::string type_name() const override;
    std::shared_ptr<SurfaceInterface> surface(int surface_id) override;
    const math::quadrature::Quadrature& integration_scheme() const override;

    // Integration-point and nodal output coordinates
    RowMatrix stress_strain_ip_rst() override;
    RowMatrix stress_strain_nodal_rst() override;

    // Quadratic triangular interpolation and MITC6-b tying
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
