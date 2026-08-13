/**
 * @file c3d8r.cpp
 * @brief Implements the reduced-integration C3D8 solid and physical hourglass stabilization.
 *
 * The one-point continuum contribution uses the common solid material-point
 * state. Hourglass control follows the Belytschko-Bindeman assumed-strain
 * stiffness form: projected Flanagan-Belytschko modes are weighted by physical
 * reference-geometry coefficients and the initial isotropic-equivalent shear
 * stiffness. No user hourglass coefficient or bulk-modulus singularity enters
 * the stabilization.
 *
 * @author Finn Eggers
 * @date 13.08.2026
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
 * Constructs the four standard Flanagan-Belytschko scalar base modes.
 *
 * The ordering `[st, tr, rs, rst]` is retained because the first three modes
 * correspond to the physical-direction indices used by the assumed-strain
 * stabilization matrix.
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
 *
 * The full two-by-two-by-two rule evaluates
 *
 *     D_bar = (1 / V0) integral D dV0,
 *
 * which defines the affine projection used by the orthogonal hourglass vectors.
 */
C3D8R::GradientMatrix C3D8R::mean_reference_gradient() {
    const auto reference_coords = node_coords_reference();
    static const math::quadrature::Quadrature full_quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_QUADRATIC
    };

    GradientMatrix integrated_gradient = GradientMatrix::Zero();
    Precision reference_volume = Precision(0);

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
        reference_volume     += measure;
    }

    logging::error(std::isfinite(reference_volume) && reference_volume > Precision(0),
        "C3D8R: invalid reference volume in element ", elem_id,
        "\nvolume: ", reference_volume);

    return integrated_gradient / reference_volume;
}

/**
 * Computes the closed-form physical geometry coefficients of the assumed-strain
 * stabilization.
 *
 * In the element-fixed reference frame the Belytschko-Bindeman coefficients are
 *
 *     H_ii = (1/3) (L_j L_k / L_i),
 *     H_ij = (1/3) L_k,                  i != j != k,
 *
 * with `L_i = Lambda_i^T x_i`. `Lambda_i` is the nodal natural-coordinate sign
 * vector and `x_i` is the corresponding nodal coordinate in the referential
 * element frame. The frame is formed from the center-Jacobian covariant axes by
 * an orientation-preserving Gram-Schmidt construction. This makes the geometry
 * coefficients invariant to a global rigid rotation while retaining the exact
 * parallelepiped expressions used by the physical stabilization.
 */
Mat3 C3D8R::hourglass_geometry_integrals() {
    const auto reference_coords = node_coords_reference();
    const auto local_coords     = node_coords_local();
    const Mat3 center_jacobian  = jacobian(
        reference_coords, Precision(0), Precision(0), Precision(0));

    Vec3 axis_1 = center_jacobian.row(0).transpose();
    Vec3 axis_2 = center_jacobian.row(1).transpose();
    Vec3 axis_3 = center_jacobian.row(2).transpose();

    logging::error(axis_1.allFinite() && axis_2.allFinite() && axis_3.allFinite(),
        "C3D8R: non-finite reference axes in element ", elem_id);
    logging::error(axis_1.norm() > Precision(0),
        "C3D8R: degenerate first reference axis in element ", elem_id);

    const Vec3 e1 = axis_1.normalized();
    axis_2 -= e1 * e1.dot(axis_2);
    logging::error(axis_2.norm() > Precision(0),
        "C3D8R: degenerate second reference axis in element ", elem_id);
    Vec3 e2 = axis_2.normalized();
    Vec3 e3 = e1.cross(e2);
    logging::error(e3.norm() > Precision(0),
        "C3D8R: degenerate third reference axis in element ", elem_id);
    e3.normalize();

    if (e3.dot(axis_3) < Precision(0)) {
        e2 = -e2;
        e3 = -e3;
    }

    StaticMatrix<D, D> frame = StaticMatrix<D, D>::Zero();
    frame.col(0) = e1;
    frame.col(1) = e2;
    frame.col(2) = e3;

    const GradientMatrix local_reference_coords = reference_coords * frame;
    StaticVector<D> generalized_lengths = StaticVector<D>::Zero();
    for (Dim i = 0; i < D; ++i) {
        generalized_lengths(i) = local_coords.col(i).dot(local_reference_coords.col(i));
        logging::error(std::isfinite(generalized_lengths(i)) && generalized_lengths(i) > Precision(0),
            "C3D8R: invalid physical hourglass length in element ", elem_id,
            "\ndirection: ", i,
            "\nlength: ", generalized_lengths(i));
    }

    Mat3 H = Mat3::Zero();
    for (Dim i = 0; i < D; ++i) {
        const Dim j = (i + 1) % D;
        const Dim k = (i + 2) % D;
        H(i, i) = Precision(1) / Precision(3) *
            generalized_lengths(j) * generalized_lengths(k) / generalized_lengths(i);
    }

    for (Dim i = 0; i < D; ++i) {
        for (Dim j = i + 1; j < D; ++j) {
            const Dim k = Precision(3) - i - j;
            H(i, j) = Precision(1) / Precision(3) * generalized_lengths(k);
            H(j, i) = H(i, j);
        }
    }

    logging::error(H.allFinite(),
        "C3D8R: non-finite physical hourglass geometry in element ", elem_id,
        "\nH: ", H);
    return H;
}

/**
 * Extracts isotropic-equivalent initial material parameters state-neutrally.
 *
 * The physical stabilization uses shear modulus `mu` and Poisson ratio `nu`.
 * They are obtained from the initial constitutive tangent as
 *
 *     mu     = mean(C44, C55, C66),
 *     lambda = mean(C12, C13, C23),
 *     nu     = lambda / (2 (lambda + mu)).
 *
 * The complete physical material-point history is restored after this auxiliary
 * zero-strain tangent query.
 */
StaticVector<2> C3D8R::hourglass_material_parameters() {
    const Index state_row = this->mp_index(0);
    Field& state_field = *this->_model_data->material_state;

    std::vector<Precision> saved_state(static_cast<std::size_t>(state_field.components));
    for (Index component = 0; component < state_field.components; ++component) {
        saved_state[static_cast<std::size_t>(component)] = state_field(state_row, component);
    }

    Precision* state = &state_field(state_row, 0);
    const Mat6 tangent = material_tangent_reference(
        Precision(0), Precision(0), Precision(0), state);

    for (Index component = 0; component < state_field.components; ++component) {
        state_field(state_row, component) = saved_state[static_cast<std::size_t>(component)];
    }

    const Precision mu =
        (tangent(3, 3) + tangent(4, 4) + tangent(5, 5)) / Precision(3);
    const Precision lame_lambda = (
        tangent(0, 1) + tangent(1, 0) +
        tangent(0, 2) + tangent(2, 0) +
        tangent(1, 2) + tangent(2, 1)
    ) / Precision(6);
    const Precision denominator = Precision(2) * (lame_lambda + mu);
    const Precision nu = lame_lambda / denominator;

    logging::error(std::isfinite(mu) && mu > Precision(0),
        "C3D8R: invalid initial shear stiffness in element ", elem_id,
        "\nmu: ", mu);
    logging::error(std::isfinite(nu) && nu > Precision(-1) && nu < Precision(0.5),
        "C3D8R: invalid isotropic-equivalent Poisson ratio in element ", elem_id,
        "\nnu: ", nu,
        "\nlambda: ", lame_lambda,
        "\nmu: ", mu);

    StaticVector<2> result = StaticVector<2>::Zero();
    result(0) = mu;
    result(1) = nu;
    return result;
}

/**
 * Builds the parameter-free Belytschko-Bindeman reference stabilization.
 *
 * The corrected vectors
 *
 *     gamma = (I - D_bar X^T) Gamma
 *
 * annihilate constant-gradient displacement fields. With `i,j,k` cyclic, the
 * eight-by-eight displacement blocks are
 *
 *     k_ii = H_ii [(1+nu)/(1-nu) (gamma_j gamma_j^T + gamma_k gamma_k^T)
 *                   + (1/3) gamma_4 gamma_4^T]
 *            + (1/2)(H_jj + H_kk) gamma_i gamma_i^T,
 *
 *     k_ij = H_ij [nu/(1-nu) gamma_j gamma_i^T
 *                   + (1/2) gamma_i gamma_j^T],
 *
 *     K_hg = 2 mu [k_ij].
 *
 * No term contains `1/(1-2 nu)`, so the stabilization remains bounded as the
 * continuum approaches incompressibility. The matrix depends only on reference
 * geometry and the initial material tangent and is cached.
 */
C3D8R::Matrix24 C3D8R::hourglass_stiffness() {
    if (hourglass_stiffness_cached) {
        return hourglass_stiffness_cache;
    }

    const auto reference_coords = node_coords_reference();
    const GradientMatrix mean_gradient = mean_reference_gradient();
    const StaticMatrix<N, N> affine_projector =
        StaticMatrix<N, N>::Identity() - mean_gradient * reference_coords.transpose();
    const HourglassModes gamma = affine_projector * primitive_hourglass_modes();

    const Mat3 H = hourglass_geometry_integrals();
    const StaticVector<2> parameters = hourglass_material_parameters();
    const Precision mu = parameters(0);
    const Precision nu = parameters(1);
    const Precision normal_factor  = (Precision(1) + nu) / (Precision(1) - nu);
    const Precision poisson_factor = nu / (Precision(1) - nu);

    Matrix24 stiffness = Matrix24::Zero();

    for (Dim i = 0; i < D; ++i) {
        const Dim j = (i + 1) % D;
        const Dim k = (i + 2) % D;

        const StaticMatrix<N, N> block =
            H(i, i) * (
                normal_factor * (
                    gamma.col(j) * gamma.col(j).transpose() +
                    gamma.col(k) * gamma.col(k).transpose()
                ) +
                Precision(1) / Precision(3) * gamma.col(3) * gamma.col(3).transpose()
            ) +
            Precision(0.5) * (H(j, j) + H(k, k)) *
                gamma.col(i) * gamma.col(i).transpose();

        for (Index node_a = 0; node_a < N; ++node_a) {
            for (Index node_b = 0; node_b < N; ++node_b) {
                stiffness(D * node_a + i, D * node_b + i) = block(node_a, node_b);
            }
        }
    }

    for (Dim i = 0; i < D; ++i) {
        for (Dim j = 0; j < D; ++j) {
            if (i == j) {
                continue;
            }

            const StaticMatrix<N, N> block = H(i, j) * (
                poisson_factor * gamma.col(j) * gamma.col(i).transpose() +
                Precision(0.5) * gamma.col(i) * gamma.col(j).transpose()
            );

            for (Index node_a = 0; node_a < N; ++node_a) {
                for (Index node_b = 0; node_b < N; ++node_b) {
                    stiffness(D * node_a + i, D * node_b + j) = block(node_a, node_b);
                }
            }
        }
    }

    hourglass_stiffness_cache = Precision(2) * mu * stiffness;
    hourglass_stiffness_cache = Precision(0.5) *
        (hourglass_stiffness_cache + hourglass_stiffness_cache.transpose());

    logging::error(hourglass_stiffness_cache.allFinite(),
        "C3D8R: non-finite physical hourglass stiffness in element ", elem_id);

    hourglass_stiffness_cached = true;
    return hourglass_stiffness_cache;
}

C3D8R::Vector24 C3D8R::local_displacement() {
    const GradientMatrix reference_coords = node_coords_reference();
    const GradientMatrix current_coords   = node_coords_current();
    const GradientMatrix displacement     = current_coords - reference_coords;

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
    mapped  = Precision(0.5) * (mapped + mapped.transpose());
    return mapped;
}

MapMatrix C3D8R::stiffness_tangent(Precision*   buffer,
                                   Field&       ip_stress_state,
                                   NodeData&    nodal_forces,
                                   const Field& displacement) {
    const Matrix24 hourglass = hourglass_stiffness();
    MapMatrix mapped = SolidElement<N>::stiffness_tangent(
        buffer, ip_stress_state, nodal_forces, displacement);

    mapped += hourglass;
    mapped  = Precision(0.5) * (mapped + mapped.transpose());
    assemble_local_force(nodal_forces, hourglass * local_displacement());
    return mapped;
}

void C3D8R::compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) {
    C3D8::compute_internal_force_nonlinear(node_forces, ip_stress);
    const Matrix24 hourglass = hourglass_stiffness();
    assemble_local_force(node_forces, hourglass * local_displacement());
}

} // namespace fem::model
