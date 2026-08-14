/**
 * @file c3d8r.cpp
 * @brief Implements the reduced-integration C3D8 solid and CalculiX-style hourglass stabilization.
 *
 * The one-point continuum contribution uses the common solid material-point
 * state. Hourglass stiffness uses the current center-point constitutive C1111
 * component state-neutrally, together with the 0.005 factor used by CalculiX.
 * The same instantaneous hourglass matrix supplies residual and tangent terms.
 *
 * @author Finn Eggers
 * @date 14.08.2026
 */

#include "c3d8r.h"

#include <cmath>
#include <vector>

namespace fem::model {

C3D8R::C3D8R(ID elem_id, const std::array<ID, N>& node_ids)
    : C3D8(elem_id, node_ids) {}

std::string C3D8R::type_name() const {
    return "C3D8R";
}

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

/**
 * Constructs the four primitive scalar hourglass modes at the C3D8 nodes.
 */
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

/**
 * Computes the volume-averaged reference shape-function gradients.
 */
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

/**
 * Evaluates the current center-point C1111 constitutive tangent state-neutrally.
 *
 * CalculiX scales C3D8R hourglass control with `stiff(1)`, i.e. the current
 * material tangent component C1111. The active material state row is saved and
 * restored around this auxiliary query so only the ordinary continuum update
 * advances constitutive history.
 */
Precision C3D8R::hourglass_material_scale() {
    const Index state_row = this->mp_index(0);
    Field& state_field = *this->_model_data->material_state;

    std::vector<Precision> saved_state(static_cast<std::size_t>(state_field.components));
    for (Index component = 0; component < state_field.components; ++component) {
        saved_state[static_cast<std::size_t>(component)] = state_field(state_row, component);
    }

    const auto reference_coords = node_coords_reference();
    const auto current_coords = node_coords_current();
    const Mat3 F = deformation_gradient(
        reference_coords, current_coords,
        Precision(0), Precision(0), Precision(0));
    const VolumeStrainGreenLagrange strain =
        VolumeStrainGreenLagrange::from_deformation_gradient(F);

    Precision* state = &state_field(state_row, 0);
    VolumeStressPK2 stress;
    Mat6 tangent;
    evaluate_material(
        Precision(0), Precision(0), Precision(0),
        strain, state, stress, tangent);

    for (Index component = 0; component < state_field.components; ++component) {
        state_field(state_row, component) = saved_state[static_cast<std::size_t>(component)];
    }

    const Precision c1111 = tangent(0, 0);
    logging::error(std::isfinite(c1111) && c1111 > Precision(0),
        "C3D8R: invalid current C1111 hourglass stiffness in element ", elem_id,
        "\nC1111: ", c1111,
        "\ndet(F): ", F.determinant());

    return c1111;
}

/**
 * Builds the instantaneous CalculiX-style hourglass stabilization tangent.
 *
 * Primitive modes are projected with `I - D_bar X^T` so affine displacement
 * fields remain unstabilized. The scalar stiffness follows the CalculiX C3D8R
 * form
 *
 *     k_hg = 0.005 C1111_current V0 sum_a ||grad_bar N_a||^2.
 *
 * The reference geometry and projected modes are unchanged; only the material
 * scale follows the current constitutive tangent. As in CalculiX, derivatives
 * of C1111 with respect to displacement are not added to the hourglass tangent.
 */
C3D8R::Matrix24 C3D8R::hourglass_stiffness() {
    const auto reference_coords = node_coords_reference();

    Precision reference_volume = Precision(0);
    const GradientMatrix mean_gradient = mean_reference_gradient(reference_volume);

    const StaticMatrix<N, N> projector =
        StaticMatrix<N, N>::Identity() - mean_gradient * reference_coords.transpose();
    const HourglassModes modes = projector * primitive_hourglass_modes();

    const Precision c1111 = hourglass_material_scale();
    const Precision gradient_scale = mean_gradient.array().square().sum();
    const Precision hourglass_scale =
        calculix_hourglass_coefficient * c1111 * reference_volume * gradient_scale;

    logging::error(std::isfinite(hourglass_scale) && hourglass_scale > Precision(0),
        "C3D8R: invalid CalculiX-style hourglass stiffness in element ", elem_id,
        "\nscale: ", hourglass_scale,
        "\nC1111: ", c1111);

    const StaticMatrix<N, N> scalar_stiffness =
        hourglass_scale * modes * modes.transpose();

    Matrix24 stiffness = Matrix24::Zero();
    for (Index node_a = 0; node_a < N; ++node_a) {
        for (Index node_b = 0; node_b < N; ++node_b) {
            for (Dim dof = 0; dof < D; ++dof) {
                stiffness(D * node_a + dof, D * node_b + dof) =
                    scalar_stiffness(node_a, node_b);
            }
        }
    }

    return Precision(0.5) * (stiffness + stiffness.transpose());
}

C3D8R::Vector24 C3D8R::local_displacement() {
    const GradientMatrix reference_coords = node_coords_reference();
    const GradientMatrix current_coords = node_coords_current();
    const GradientMatrix displacement = current_coords - reference_coords;

    Vector24 result = Vector24::Zero();
    for (Index node = 0; node < N; ++node) {
        for (Dim dof = 0; dof < D; ++dof) {
            result(D * node + dof) = displacement(node, dof);
        }
    }
    return result;
}

void C3D8R::assemble_local_force(Field& node_forces, const Vector24& local_force) {
    logging::error(node_forces.domain == FieldDomain::NODE,
        "C3D8R: internal force output must use NODE domain");
    logging::error(node_forces.components >= D,
        "C3D8R: internal force output requires at least three components");

    for (Index node = 0; node < N; ++node) {
        const Index node_id = static_cast<Index>(node_ids[node]);
        for (Dim dof = 0; dof < D; ++dof) {
            node_forces(node_id, dof) += local_force(D * node + dof);
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

void C3D8R::compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) {
    C3D8::compute_internal_force_nonlinear(node_forces, ip_stress);

    const Matrix24 hourglass = hourglass_stiffness();
    assemble_local_force(node_forces, hourglass * local_displacement());
}

} // namespace fem::model
