/**
 * @file c3d8r.cpp
 * @brief Implements the reduced-integration C3D8 solid and hourglass stabilization.
 *
 * Hourglass scaling queries the state-neutral elastic backbone rather than the
 * trial constitutive law, so plastic history is never touched by stabilization.
 *
 * @author Finn Eggers
 * @date 16.08.2026
 */

#include "c3d8r.h"

#include <cmath>

namespace fem::model {

C3D8R::C3D8R(ID elem_id, const std::array<ID, N>& node_ids)
    : C3D8(elem_id, node_ids) {}

std::string C3D8R::type_name() const { return "C3D8R"; }

const math::quadrature::Quadrature& C3D8R::integration_scheme_stiffness() const {
    static const math::quadrature::Quadrature quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_CONSTANT
    };
    return quadrature;
}

RowMatrix C3D8R::stress_strain_nodal_rst() {
    return RowMatrix::Zero(N, D);
}

C3D8R::HourglassModes C3D8R::primitive_hourglass_modes() {
    const auto local_coords = node_coords_local();
    HourglassModes modes = HourglassModes::Zero();
    for (Index node = 0; node < N; ++node) {
        const Precision r = local_coords(node, 0);
        const Precision s = local_coords(node, 1);
        const Precision t = local_coords(node, 2);
        modes(node, 0) = s * t;
        modes(node, 1) = t * r;
        modes(node, 2) = r * s;
        modes(node, 3) = r * s * t;
    }
    return modes;
}

C3D8R::GradientMatrix C3D8R::mean_reference_gradient(Precision& reference_volume) {
    const auto reference_coords = node_coords_reference();
    static const math::quadrature::Quadrature full_quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_QUADRATIC
    };

    GradientMatrix integrated_gradient = GradientMatrix::Zero();
    reference_volume = Precision(0);
    for (Index q = 0; q < full_quadrature.count(); ++q) {
        const auto point = full_quadrature.get_point(q);
        Precision det0 = Precision(0);
        const auto gradient = shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        logging::error(std::isfinite(det0) && det0 > Precision(0),
            "C3D8R: invalid reference determinant in element ", elem_id,
            "\ndet(J0): ", det0);
        const Precision measure = det0 * point.w;
        integrated_gradient += gradient * measure;
        reference_volume += measure;
    }

    logging::error(std::isfinite(reference_volume) && reference_volume > Precision(0),
        "C3D8R: invalid reference volume in element ", elem_id,
        "\nvolume: ", reference_volume);
    return integrated_gradient / reference_volume;
}

Precision C3D8R::hourglass_material_scale() {
    const Mat6 material_tangent = material_tangent_reference(
        Precision(0), Precision(0), Precision(0));
    const Precision shear_scale =
        (material_tangent(3, 3)
       + material_tangent(4, 4)
       + material_tangent(5, 5)) / Precision(3);
    logging::error(std::isfinite(shear_scale) && shear_scale > Precision(0),
        "C3D8R: invalid initial mean shear stiffness in element ", elem_id,
        "\nscale: ", shear_scale);
    return shear_scale;
}

C3D8R::Matrix24 C3D8R::hourglass_stiffness() {
    const auto reference_coords = node_coords_reference();
    Precision reference_volume = Precision(0);
    const GradientMatrix mean_gradient = mean_reference_gradient(reference_volume);
    const StaticMatrix<N, N> projector =
        StaticMatrix<N, N>::Identity() - mean_gradient * reference_coords.transpose();
    const HourglassModes modes = projector * primitive_hourglass_modes();

    const Precision material_scale = hourglass_material_scale();
    const Precision gradient_scale = mean_gradient.array().square().sum();
    const Precision hourglass_scale = default_hourglass_coefficient
                                    * material_scale
                                    * reference_volume
                                    * gradient_scale;
    logging::error(std::isfinite(hourglass_scale) && hourglass_scale > Precision(0),
        "C3D8R: invalid hourglass stiffness in element ", elem_id,
        "\nscale: ", hourglass_scale);

    const StaticMatrix<N, N> scalar_stiffness =
        hourglass_scale * modes * modes.transpose();
    Matrix24 stiffness = Matrix24::Zero();
    for (Index a = 0; a < N; ++a) {
        for (Index b = 0; b < N; ++b) {
            for (Dim d = 0; d < D; ++d) {
                stiffness(D * a + d, D * b + d) = scalar_stiffness(a, b);
            }
        }
    }
    return Precision(0.5) * (stiffness + stiffness.transpose());
}

C3D8R::Vector24 C3D8R::local_displacement() {
    const GradientMatrix displacement =
        node_coords_current() - node_coords_reference();
    Vector24 result = Vector24::Zero();
    for (Index node = 0; node < N; ++node) {
        for (Dim d = 0; d < D; ++d) {
            result(D * node + d) = displacement(node, d);
        }
    }
    return result;
}

void C3D8R::assemble_local_force(Field& node_forces,
                                 const Vector24& local_force) {
    logging::error(node_forces.domain == FieldDomain::NODE,
        "C3D8R: internal force output must use NODE domain");
    logging::error(node_forces.components >= D,
        "C3D8R: internal force output requires at least three components");
    for (Index node = 0; node < N; ++node) {
        const Index node_id = static_cast<Index>(node_ids[node]);
        for (Dim d = 0; d < D; ++d) {
            node_forces(node_id, d) += local_force(D * node + d);
        }
    }
}

MapMatrix C3D8R::stiffness(Precision* buffer) {
    MapMatrix mapped{buffer, ndof, ndof};
    C3D8::stiffness(buffer);
    mapped += hourglass_stiffness();
    mapped = Precision(0.5) * (mapped + mapped.transpose());
    return mapped;
}

MapMatrix C3D8R::stiffness_tangent(Precision* buffer,
                                   Field& ip_stress_state,
                                   NodeData& nodal_forces,
                                   const Field& displacement) {
    const Matrix24 hourglass = hourglass_stiffness();
    MapMatrix mapped = SolidElement<N>::stiffness_tangent(
        buffer, ip_stress_state, nodal_forces, displacement);
    mapped += hourglass;
    mapped = Precision(0.5) * (mapped + mapped.transpose());
    assemble_local_force(nodal_forces, hourglass * local_displacement());
    return mapped;
}

void C3D8R::compute_internal_force_nonlinear(Field& node_forces,
                                             const Field& ip_stress) {
    C3D8::compute_internal_force_nonlinear(node_forces, ip_stress);
    const Matrix24 hourglass = hourglass_stiffness();
    assemble_local_force(node_forces, hourglass * local_displacement());
}

} // namespace fem::model
