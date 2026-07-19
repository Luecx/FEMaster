/**
 * @file c3d8r.cpp
 * @brief Implements the reduced-integration C3D8 solid with hourglass control.
 */

#include "c3d8r.h"

#include <algorithm>
#include <cmath>

namespace fem::model {

C3D8R::C3D8R(ID elem_id, const std::array<ID, 8>& node_ids)
    : C3D8(elem_id, node_ids) {}

const math::quadrature::Quadrature& C3D8R::integration_scheme() const {
    static const math::quadrature::Quadrature quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_CONSTANT
    };
    return quadrature;
}

RowMatrix C3D8R::stress_strain_nodal_rst() {
    return RowMatrix::Zero(8, 3);
}

C3D8R::Matrix24 C3D8R::hourglass_stiffness() {
    // Compute volume-averaged reference gradients
    const auto reference_coords = node_coords_reference();

    static const math::quadrature::Quadrature full_quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_QUADRATIC
    };

    StaticMatrix<8, 3> mean_gradient    = StaticMatrix<8, 3>::Zero();
    Precision          reference_volume = Precision(0);

    for (Index q = 0; q < full_quadrature.count(); ++q) {
        const auto point = full_quadrature.get_point(q);

        Precision det0 = Precision(0);
        const auto gradient = shape_derivatives_reference(
            reference_coords,
            point.r,
            point.s,
            point.t,
            det0
        );
        const Precision measure = det0 * point.w;

        mean_gradient     += gradient * measure;
        reference_volume += measure;
    }

    logging::error(std::isfinite(reference_volume) && reference_volume > Precision(0),
                   "C3D8R: invalid reference volume in element ", elem_id);

    mean_gradient /= reference_volume;

    // Build and project the four primitive hourglass modes
    const auto local_coords = node_coords_local();
    HourglassModes primitive = HourglassModes::Zero();

    for (Index node = 0; node < 8; ++node) {
        const Precision r = local_coords(node, 0);
        const Precision s = local_coords(node, 1);
        const Precision t = local_coords(node, 2);

        primitive(node, 0) = s * t;
        primitive(node, 1) = t * r;
        primitive(node, 2) = r * s;
        primitive(node, 3) = r * s * t;
    }

    const StaticMatrix<8, 8> projector =
        StaticMatrix<8, 8>::Identity() - mean_gradient * reference_coords.transpose();
    const HourglassModes modes = projector * primitive;

    // Scale the stabilization with the largest normal material stiffness
    const Mat6 material_tangent = material_tangent_reference(
        Precision(0),
        Precision(0),
        Precision(0)
    );
    const Precision material_scale  = std::max({
        material_tangent(0, 0),
        material_tangent(1, 1),
        material_tangent(2, 2)
    });
    const Precision gradient_scale  = mean_gradient.array().square().sum();
    const Precision hourglass_scale = Precision(0.005) * material_scale * reference_volume * gradient_scale;

    logging::error(std::isfinite(material_scale) && material_scale > Precision(0),
                   "C3D8R: invalid material stiffness scale in element ", elem_id);
    logging::error(std::isfinite(hourglass_scale) && hourglass_scale > Precision(0),
                   "C3D8R: invalid hourglass stiffness in element ", elem_id);

    const StaticMatrix<8, 8> scalar_stiffness = hourglass_scale * modes * modes.transpose();

    Matrix24 hourglass_stiffness = Matrix24::Zero();
    for (Index node_a = 0; node_a < 8; ++node_a) {
        for (Index node_b = 0; node_b < 8; ++node_b) {
            for (Dim dof = 0; dof < 3; ++dof) {
                hourglass_stiffness(3 * node_a + dof, 3 * node_b + dof) =
                    scalar_stiffness(node_a, node_b);
            }
        }
    }

    return hourglass_stiffness;
}

MapMatrix C3D8R::stiffness(Precision* buffer) {
    MapMatrix stiffness{buffer, 24, 24};

    C3D8::stiffness(buffer);
    stiffness += hourglass_stiffness();

    return stiffness;
}

void C3D8R::compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) {
    C3D8::compute_internal_force_nonlinear(node_forces, ip_stress);

    // Add the Newton-consistent reference hourglass force
    const auto reference_coords = node_coords_reference();
    const auto current_coords   = node_coords_current();

    Vector24 displacement = Vector24::Zero();
    for (Index node = 0; node < 8; ++node) {
        for (Dim dof = 0; dof < 3; ++dof) {
            displacement(3 * node + dof) = current_coords(node, dof) - reference_coords(node, dof);
        }
    }

    const Vector24 hourglass_force = hourglass_stiffness() * displacement;

    for (Index node = 0; node < 8; ++node) {
        const Index node_id = static_cast<Index>(node_ids[node]);
        for (Dim dof = 0; dof < 3; ++dof) {
            node_forces(node_id, dof) += hourglass_force(3 * node + dof);
        }
    }
}

} // namespace fem::model
