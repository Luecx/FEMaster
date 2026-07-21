/**
 * @file frt_shell_s8.h
 * @brief Declares the eight-node finite-rotation MITC8 shell element.
 *
 * The element combines the common Total-Lagrangian finite-rotation shell
 * formulation with the classical MITC8 assumed in-layer and transverse-shear
 * strain fields. The quadratic serendipity geometry is evaluated directly on
 * the curved reference midsurface.
 *
 * @see FRTShell
 * @see Surface8
 *
 * @author Finn Eggers
 * @date 20.07.2026
 */

#pragma once

#include "frt_shell.h"
#include "../geometry/surface/surface8.h"

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace fem::model {

/**
 * @brief Eight-node quadratic finite-rotation MITC8 shell element.
 *
 * The first four nodes are the quadrilateral corners and the final four nodes
 * are the midside nodes on edges `0-1`, `1-2`, `2-3` and `3-0`. The element
 * uses the classical MITC8 interpolation for both in-layer and transverse-
 * shear strains. Membrane strains and curvatures use the same in-layer tensor
 * interpolation, while the two transverse-shear components use separate
 * five-function assumed fields.
 */
struct FRTShellS8 : FRTShell<8> {
    // Surface geometry and numerical integration
    Surface8                     geometry;
    math::quadrature::Quadrature integration_scheme_;

    // Construction
    FRTShellS8(ID id, const std::array<ID, 8>& nodes);
    ~FRTShellS8() override = default;

    // Element identification and surface connectivity
    std::string type_name() const override;
    std::shared_ptr<SurfaceInterface> surface(int surface_id) override;

    // Full 3 x 3 integration on the quadratic reference surface
    const math::quadrature::Quadrature& integration_scheme() const override;
    RowMatrix stress_strain_ip_rst   () override;
    RowMatrix stress_strain_nodal_rst() override;

    // Quadratic serendipity interpolation and natural node coordinates
    VecN  shape_function     (Precision r, Precision s) const override;
    MatN2 shape_derivative   (Precision r, Precision s) const override;
    MatN2 node_coords_natural() const override;

    // Classical MITC8 in-layer and transverse-shear tying points. The assumed
    // tensor fields are reconstructed in physical reference-surface bases and
    // finally returned in natural covariant components.
    std::vector<Vec2> tying_point_coordinates() const override;

    void apply_mitc_natural(
        const EvaluationData& data,
        const ReferencePoint& point,
        Vec8&                 strain_nat,
        Mat8x6N*              B_nat
    ) const override;

    void pull_back_mitc_resultants(
        const ReferencePoint& point,
        const Vec8&           assumed_weights,
        Vec8&                 compatible_weights,
        std::vector<Vec8>&    tying_weights
    ) const override;
};

} // namespace fem::model
