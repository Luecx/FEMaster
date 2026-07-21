/**
 * @file frt_shell_s3.h
 * @brief Declares the three-node finite-rotation MITC3 shell element.
 *
 * The element combines the common Total-Lagrangian finite-rotation shell
 * kernel with the isotropic MITC3 transverse-shear interpolation. The geometry
 * is a linear triangular midsurface, and the three tying points lie at the
 * edge midpoints of the natural triangle.
 *
 * The implementation uses the original MITC3 assumed shear field. It does not
 * introduce the internal cubic bubble rotations of MITC3+.
 *
 * @see FRTShell
 * @see Surface3
 *
 * @author Finn Eggers
 * @date 20.07.2026
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
    // Element geometry and triangular area quadrature
    Surface3                     geometry;
    math::quadrature::Quadrature integration_scheme_;

    // Construction
    FRTShellS3(ID id, const std::array<ID, 3>& nodes);
    ~FRTShellS3() override = default;

    // Element identification, surface extraction and numerical integration
    std::string type_name() const override;
    std::shared_ptr<SurfaceInterface> surface(int surface_id) override;
    const math::quadrature::Quadrature& integration_scheme() const override;

    // Integration-point and nodal output coordinates
    RowMatrix stress_strain_ip_rst() override;
    RowMatrix stress_strain_nodal_rst() override;

    // Linear triangular interpolation and MITC3 tying
    VecN  shape_function     (Precision r, Precision s) const override;
    MatN2 shape_derivative   (Precision r, Precision s) const override;
    MatN2 node_coords_natural() const override;
    std::vector<Vec2> tying_point_coordinates() const override;

    void apply_mitc_natural(const EvaluationData& data,
                            const ReferencePoint& point,
                            Vec8&                 strain_nat,
                            Mat8x6N*              B_nat) const override;

    void pull_back_mitc_resultants(const ReferencePoint& point,
                                   const Vec8&           assumed_weights,
                                   Vec8&                 compatible_weights,
                                   std::vector<Vec8>&    tying_weights) const override;
};

} // namespace fem::model
