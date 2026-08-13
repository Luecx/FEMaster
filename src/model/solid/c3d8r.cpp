/**
 * @file c3d8r.cpp
 * @brief Implements the reduced-integration C3D8 solid and physical hourglass stabilization.
 *
 * The one-point continuum contribution uses the common solid material-point
 * state. Hourglass control follows the Belytschko-Bindeman assumed-strain
 * stiffness form: projected Flanagan-Belytschko modes are weighted by exact
 * reference-geometry integrals and the initial isotropic-equivalent shear
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
    // One-point hexahedral integration at r = s = t = 0 with weight eight.
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
 *
 * The natural-coordinate products `[s t, t r, r s, r s t]` are the four
 * trilinear fields omitted by one-point integration. This ordering is the
 * physical-stabilization convention `[h1, h2, h3, h4]`, where `h1` omits `r`,
 * `h2` omits `s`, and `h3` omits `t`.
 *
 * @return Eight-by-four primitive modal matrix in element-node ordering.
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
 * A full two-by-two-by-two rule evaluates
 *
 *     D_bar = (1 / V0) integral D dV0.
 *
 * These mean gradients define the Flanagan-Belytschko correction that removes
 * all affine displacement fields from the primitive hourglass vectors.
 *
 * @return Mean reference gradient with one row per element node.
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
            reference_coords,
            point.r,
            point.s,
            point.t,
            det0
        );

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
 * Computes the geometry coefficients of the assumed-strain stabilization.
 *
 * The natural hourglass fields are
 *
 *     h1 = s t,  h2 = t r,  h3 = r s,  h4 = r s t.
 *
 * Belytschko-Bindeman physical stabilization uses the length-valued coefficients
 *
 *     H_ii = integral (h_j,i)^2 dV
 *          = integral (h_k,i)^2 dV
 *          = 3 integral (h_4,i)^2 dV,
 *
 *     H_ij = integral h_i,j h_j,i dV.
 *
 * The three equivalent expressions for `H_ii` coincide for the associated
 * parallelepiped. For a general trilinear C3D8 geometry their symmetric average
 * is integrated numerically with the full two-by-two-by-two rule; this avoids a
 * coordinate-direction preference while retaining the physical dimensions and
 * exact parallelepiped limit of the assumed-strain formula.
 *
 * @return Symmetric three-by-three matrix of physical geometry coefficients.
 */
Mat3 C3D8R::hourglass_geometry_integrals() {
    const auto reference_coords = node_coords_reference();

    static const math::quadrature::Quadrature full_quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_QUADRATIC
    };

    Mat3 H = Mat3::Zero();

    for (Index q = 0; q < full_quadrature.count(); ++q) {
        const auto point = full_quadrature.get_point(q);
        const Mat3 J0 = jacobian(reference_coords, point.r, point.s, point.t);
        const Precision det0 = J0.determinant();

        logging::error(std::isfinite(det0) && det0 > Precision(0),
            "C3D8R: invalid reference determinant in hourglass geometry for element ", elem_id,
            "\ndet(J0): ", det0);

        StaticMatrix<4, D> dh_dxi = StaticMatrix<4, D>::Zero();
        dh_dxi(0, 1) = point.t;
        dh_dxi(0, 2) = point.s;
        dh_dxi(1, 0) = point.t;
        dh_dxi(1, 2) = point.r;
        dh_dxi(2, 0) = point.s;
        dh_dxi(2, 1) = point.r;
        dh_dxi(3, 0) = point.s * point.t;
        dh_dxi(3, 1) = point.r * point.t;
        dh_dxi(3, 2) = point.r * point.s;

        const StaticMatrix<4, D> dh_dx =
            (J0.inverse() * dh_dxi.transpose()).transpose();
        const Precision measure = det0 * point.w;

        for (Dim i = 0; i < D; ++i) {
            const Dim j = (i + 1) % D;
            const Dim k = (i + 2) % D;

            const Precision hii = (
                dh_dx(j, i) * dh_dx(j, i) +
                dh_dx(k, i) * dh_dx(k, i) +
                Precision(3) * dh_dx(3, i) * dh_dx(3, i)
            ) / Precision(3);

            H(i, i) += hii * measure;

            for (Dim m = i + 1; m < D; ++m) {
                const Precision hij = dh_dx(i, m) * dh_dx(m, i) * measure;
                H(i, m) += hij;
                H(m, i) += hij;
            }
        }
    }

    logging::error(H.allFinite(),
        "C3D8R: non-finite physical hourglass geometry in element ", elem_id);
    logging::error(H(0, 0) > Precision(0) && H(1, 1) > Precision(0) && H(2, 2) > Precision(0),
        "C3D8R: non-positive physical hourglass geometry in element ", elem_id,
        "\nH: ", H);

    return H;
}

/**
 * Extracts isotropic-equivalent initial material parameters state-neutrally.
 *
 * The Belytschko-Bindeman solid stabilization is written in terms of shear
 * modulus `mu` and Poisson ratio `nu`. FEMaster obtains both from the initial
 * six-by-six constitutive tangent so hyperelastic materials use their exact
 * infinitesimal tangent without a separate material-type dependency:
 *
 *     mu     = mean(C44, C55, C66),
 *     lambda = mean(C12, C13, C23),
 *     nu     = lambda / (2 (lambda + mu)).
 *
 * The complete material-point state row is restored after the auxiliary tangent
 * evaluation. For an isotropic or initially isotropic hyperelastic law these are
 * the exact initial Lamé-equivalent parameters; for a more general tangent they
 * are its isotropic projection used only by the stabilization.
 *
 * @return Vector `[mu, nu]`.
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
 * The corrected scalar modes are
 *
 *     gamma = (I - D_bar X^T) Gamma,
 *
 * so constant-gradient displacement fields remain exactly unstabilized. With
 * `i,j,k` denoting cyclic physical directions, the eight-by-eight displacement
 * blocks follow the assumed-strain physical stabilization form
 *
 *     k_ii = H_ii [ (1+nu)/(1-nu) (gamma_j gamma_j^T + gamma_k gamma_k^T)
 *                    + 1/3 gamma_4 gamma_4^T ]
 *            + 1/2 (H_jj + H_kk) gamma_i gamma_i^T,
 *
 *     k_ij = H_ij [ nu/(1-nu) gamma_j gamma_i^T
 *                    + 1/2 gamma_i gamma_j^T ],
 *
 *     K_hg = 2 mu [k_ij].
 *
 * There is no factor containing `1/(1-2 nu)`, so the stabilization remains
 * bounded as the continuum material approaches incompressibility. No empirical
 * hourglass coefficient is used. The matrix depends only on reference geometry
 * and the initial material tangent and is therefore cached.
 *
 * @return Constant 24-by-24 physical hourglass tangent.
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
    const Precision one_minus_nu = Precision(1) - nu;
    const Precision normal_factor = (Precision(1) + nu) / one_minus_nu;
    const Precision poisson_factor = nu / one_minus_nu;

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
                Precision(1) / Precision(3) *
                    gamma.col(3) * gamma.col(3).transpose()
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

/**
 * Collects current element translations relative to the reference geometry.
 *
 * @return Twenty-four-component displacement vector in node-major XYZ ordering.
 */
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

/**
 * Scatters one element-local translational force vector into the global nodal
 * force field.
 *
 * @param node_forces Global nodal accumulator with at least XYZ components.
 * @param local_force Element force in node-major XYZ ordering.
 */
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

/**
 * Assembles the one-point material stiffness with physical hourglass control.
 *
 * @param buffer Caller-provided dense 24-by-24 storage.
 * @return Mapped symmetric continuum-plus-hourglass stiffness.
 */
MapMatrix C3D8R::stiffness(Precision* buffer) {
    MapMatrix mapped{buffer, ndof, ndof};

    C3D8::stiffness(buffer);
    mapped += hourglass_stiffness();
    mapped  = Precision(0.5) * (mapped + mapped.transpose());

    return mapped;
}

/**
 * Assembles the one-point continuum tangent and matching physical hourglass
 * contribution.
 *
 * The stabilization matrix is reference-constant. Consequently the force
 * `K_hg u_e` and tangent `K_hg` are exact derivatives of the same quadratic
 * stabilization energy and do not introduce a missing material-derivative term.
 *
 * @param buffer Caller-provided dense 24-by-24 tangent storage.
 * @param ip_stress_state Global center-point PK2 stress field to update.
 * @param nodal_forces Global nodal internal-force field to increment.
 * @param displacement Trial displacement defining the current configuration.
 * @return Mapped consistent continuum-plus-hourglass tangent.
 */
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

/**
 * Recovers nonlinear internal force from stored continuum stress and the
 * current physical hourglass displacement.
 *
 * @param node_forces Global nodal internal-force field to increment.
 * @param ip_stress Stored center-point second Piola-Kirchhoff stress field.
 */
void C3D8R::compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) {
    C3D8::compute_internal_force_nonlinear(node_forces, ip_stress);

    const Matrix24 hourglass = hourglass_stiffness();
    assemble_local_force(node_forces, hourglass * local_displacement());
}

} // namespace fem::model
