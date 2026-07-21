/**
 * @file frt_shell_output.inl
 * @brief Implements generalized and physical shell result recovery.
 *
 * The routines evaluate MITC generalized strains and shell resultants at
 * arbitrary natural coordinates. Physical through-thickness strain and stress
 * tensors are reconstructed from membrane resultants, bending moments and
 * transverse shear resultants. Nonlinear stress output transforms the recovered
 * second Piola-Kirchhoff stress into Cauchy stress through the shell deformation
 * gradient.
 *
 * @see FRTShell
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#include "frt_shell.h"

#include "../../core/logging.h"
#include "../../material/strain/volume_strain_green_lagrange.h"
#include "../../material/stress/volume_stress_cauchy.h"
#include "../../material/stress/volume_stress_pk2.h"
#include "../../math/vec_util.h"

#include <cmath>

namespace fem::model {

using math::normalized;

/**
 * Evaluates the generalized shell strain at an arbitrary natural point.
 *
 * In nonlinear mode the strain is evaluated directly from the supplied current
 * state. In linear mode the reference B matrix is evaluated at the point and
 * multiplied by the element displacement vector.
 *
 * @param data Active evaluation data containing state and MITC tying values.
 * @param q Element displacement vector used by linear recovery.
 * @param r First natural output coordinate.
 * @param s Second natural output coordinate.
 * @param nonlinear Select direct nonlinear or linearized strain recovery.
 * @return Generalized local shell strain vector.
 */
template<Index N>
typename FRTShell<N>::Vec8 FRTShell<N>::generalized_strain_at(
    const EvaluationData& data,
    const Vec6N&          q,
    Precision             r,
    Precision             s,
    bool                  nonlinear
) const {
    const ReferencePoint* cached = cached_reference_point(r, s);
    const ReferencePoint  temporary = cached ? ReferencePoint{} : make_reference_point(r, s, Precision(0));
    const ReferencePoint& point = cached ? *cached : temporary;

    Vec8     strain_nat;
    Mat8x6N B_nat;

    if (nonlinear) {
        compute_natural_strain(data, point, strain_nat);
        apply_mitc_natural(data, point, strain_nat, nullptr);
        transform_strain_to_local(point, strain_nat);
        return strain_nat;
    }

    compute_natural_strain(data, point, strain_nat, &B_nat);
    apply_mitc_natural(data, point, strain_nat, &B_nat);
    transform_strain_to_local(point, strain_nat, &B_nat);
    return B_nat * q;
}

/**
 * Evaluates generalized shell resultants at an arbitrary natural point.
 *
 * Nonlinear output calls the shell section with the recovered generalized
 * strain. Linear output multiplies the zero-strain local section tangent by the
 * linearized generalized strain.
 *
 * @param data Active evaluation data containing state and MITC tying values.
 * @param q Element displacement vector used by linear recovery.
 * @param r First natural output coordinate.
 * @param s Second natural output coordinate.
 * @param nonlinear Select nonlinear section evaluation or linear multiplication.
 * @param strain_out Optional generalized strain output.
 * @return Generalized local shell resultant vector.
 */
template<Index N>
typename FRTShell<N>::Vec8 FRTShell<N>::generalized_resultant_at(
    const EvaluationData& data,
    const Vec6N&          q,
    Precision             r,
    Precision             s,
    bool                  nonlinear,
    Vec8*                 strain_out
) const {
    const Vec8 strain_values = generalized_strain_at(data, q, r, s, nonlinear);

    if (strain_out) {
        *strain_out = strain_values;
    }

    if (!nonlinear) {
        return resultant_stiffness(r, s) * strain_values;
    }

    ShellGeneralizedStrain strain(strain_values);
    ShellStressResultants  resultants;
    Mat8                   tangent;

    shell_section()->evaluate(
        reference_position(r, s),
        reference_basis_global(r, s),
        strain,
        true,
        resultants,
        tangent
    );

    return topology_stiffness_scale() * resultants.values();
}

/**
 * Evaluates the shell deformation gradient at one through-thickness point.
 *
 * The reference and current covariant bases include the linear director
 * variation through the thickness coordinate `z`. The result maps the
 * undeformed shell basis into the current shell basis.
 *
 * @param state Current nodal shell state.
 * @param r First natural output coordinate.
 * @param s Second natural output coordinate.
 * @param z Physical through-thickness coordinate measured from the midsurface.
 * @return Three-dimensional shell deformation gradient.
 */
template<Index N>
Mat3 FRTShell<N>::deformation_gradient_at(const CurrentState& state,
                                          Precision           r,
                                          Precision           s,
                                          Precision           z) const {
    const ReferencePoint point = make_reference_point(r, s, Precision(0));

    Vec3 x_a = Vec3::Zero();
    Vec3 x_b = Vec3::Zero();
    Vec3 d   = Vec3::Zero();
    Vec3 d_a = Vec3::Zero();
    Vec3 d_b = Vec3::Zero();

    for (Index node = 0; node < num_nodes; ++node) {
        const Vec3 x_i = state.x.row(node).transpose();
        const Vec3 d_i = state.d.row(node).transpose();

        x_a += point.dshape_ab.col(0)(node) * x_i;
        x_b += point.dshape_ab.col(1)(node) * x_i;
        d   += point.shape(node)    * d_i;
        d_a += point.dshape_ab.col(0)(node) * d_i;
        d_b += point.dshape_ab.col(1)(node) * d_i;
    }

    Mat3 reference_covariant;
    reference_covariant.col(0) = point.X_ab.col(0) + z * point.D_a;
    reference_covariant.col(1) = point.X_ab.col(1) + z * point.D_b;
    reference_covariant.col(2) = point.D;

    Mat3 current_covariant;
    current_covariant.col(0) = x_a + z * d_a;
    current_covariant.col(1) = x_b + z * d_b;
    current_covariant.col(2) = d;

    const Precision reference_det = reference_covariant.determinant();
    logging::error(std::abs(reference_det) > Precision(1e-14),
                   "FRTShell: singular reference shell basis during stress recovery in element ",
                   this->elem_id);

    return current_covariant * reference_covariant.inverse();
}

/**
 * Reconstructs physical through-thickness strain and stress tensors.
 *
 * Membrane strain and curvature are combined linearly through the thickness.
 * Membrane resultants and moments are converted into the corresponding stress
 * distribution for a homogeneous shell section. Transverse shear is assumed
 * constant through the thickness, consistent with the Reissner-Mindlin model.
 *
 * @param data Active evaluation data and current nodal state.
 * @param q Element displacement vector used by linear recovery.
 * @param r First natural output coordinate.
 * @param s Second natural output coordinate.
 * @param zeta Normalized through-thickness coordinate in `[-1,1]`.
 * @param nonlinear Select nonlinear PK2-to-Cauchy transformation.
 * @param strain_out Reconstructed global Green-Lagrange strain vector.
 * @param stress_out Reconstructed global Cauchy stress vector.
 */
template<Index N>
void FRTShell<N>::physical_stress_strain_at(
    const EvaluationData& data,
    const Vec6N&          q,
    Precision             r,
    Precision             s,
    Precision             zeta,
    bool                  nonlinear,
    Vec6&                 strain_out,
    Vec6&                 stress_out
) const {
    Vec8 generalized_strain;
    const Vec8 resultants = generalized_resultant_at(
        data,
        q,
        r,
        s,
        nonlinear,
        &generalized_strain
    );

    const Precision h = this->get_section()->thickness_;
    const Precision z = Precision(0.5) * h * zeta;

    const Index membrane_strain_start = static_cast<Index>(ShellGeneralizedStrain::Component::EpsilonXX);
    const Index curvature_start       = static_cast<Index>(ShellGeneralizedStrain::Component::KappaXX);
    const Index shear_strain_start    = static_cast<Index>(ShellGeneralizedStrain::Component::GammaXZ);
    const Index membrane_stress_start = static_cast<Index>(ShellStressResultants::Component::NXX);
    const Index moment_start          = static_cast<Index>(ShellStressResultants::Component::MXX);
    const Index shear_stress_start    = static_cast<Index>(ShellStressResultants::Component::QX);

    const Vec3 plane_strain =
        generalized_strain.template segment<3>(membrane_strain_start)
        + z * generalized_strain.template segment<3>(curvature_start);

    const Vec3 plane_stress =
        resultants.template segment<3>(membrane_stress_start) / h
        + z * (Precision(12) / (h * h * h))
        * resultants.template segment<3>(moment_start);

    const Vec2 shear_strain = generalized_strain.template segment<2>(shear_strain_start);
    const Vec2 shear_stress = resultants.template segment<2>(shear_stress_start) / h;

    VolumeStrainGreenLagrange strain_local;
    VolumeStressPK2           stress_local;

    strain_local[VolumeStrain::Component::XX]      = plane_strain(0);
    strain_local[VolumeStrain::Component::YY]      = plane_strain(1);
    strain_local[VolumeStrain::Component::GammaYZ] = shear_strain(1);
    strain_local[VolumeStrain::Component::GammaXZ] = shear_strain(0);
    strain_local[VolumeStrain::Component::GammaXY] = plane_strain(2);

    stress_local[VolumeStress::Component::XX] = plane_stress(0);
    stress_local[VolumeStress::Component::YY] = plane_stress(1);
    stress_local[VolumeStress::Component::YZ] = shear_stress(1);
    stress_local[VolumeStress::Component::XZ] = shear_stress(0);
    stress_local[VolumeStress::Component::XY] = plane_stress(2);

    const Mat3 reference_basis = reference_basis_global(r, s);
    const Mat3 green_lagrange_global = reference_basis
                                     * strain_local.tensor()
                                     * reference_basis.transpose();
    const Mat3 second_pk_global = reference_basis
                                * stress_local.tensor()
                                * reference_basis.transpose();

    strain_out = VolumeStrainGreenLagrange(green_lagrange_global).voigt();

    if (!nonlinear) {
        stress_out = VolumeStressCauchy(second_pk_global).voigt();
        return;
    }

    const Mat3      F = deformation_gradient_at(data.state, r, s, z);
    const Precision J = F.determinant();

    logging::error(J > Precision(0) && std::isfinite(J),
                   "FRTShell: invalid deformation gradient during stress recovery in element ",
                   this->elem_id,
                   ", J = ", J);

    const Mat3 cauchy_global = (F * second_pk_global * F.transpose()) / J;
    stress_out = VolumeStressCauchy(cauchy_global).voigt();
}

/**
 * Computes physical stress and strain tensors at requested natural and
 * through-thickness coordinates.
 *
 * @param strain Optional physical strain output field.
 * @param stress Optional physical stress output field.
 * @param displacement Global nodal displacement field.
 * @param rst Requested natural and normalized thickness coordinates.
 * @param offset First output row belonging to this element.
 * @param use_green_lagrange_nl Select nonlinear or linearized recovery.
 */
template<Index N>
void FRTShell<N>::compute_stress_strain(Field*           strain,
                                        Field*           stress,
                                        const Field&     displacement,
                                        const RowMatrix& rst,
                                        int              offset,
                                        bool             use_green_lagrange_nl) {
    logging::error(strain != nullptr || stress != nullptr,
                   "FRTShell: compute_stress_strain requires at least one output field");
    logging::error(rst.cols() >= 3,
                   "FRTShell: stress/strain coordinates require r, s and t columns");

    const CurrentState state = use_green_lagrange_nl
        ? current_state_from_displacement(displacement)
        : reference_state();

    const EvaluationData data = init_evaluation(
        state,
        true,
        !use_green_lagrange_nl,
        false,
        false
    );

    const Vec6N q = element_displacement_vector(displacement);

    for (Eigen::Index point = 0; point < rst.rows(); ++point) {
        Vec6 strain_value;
        Vec6 stress_value;

        physical_stress_strain_at(
            data,
            q,
            rst(point, 0),
            rst(point, 1),
            rst(point, 2),
            use_green_lagrange_nl,
            strain_value,
            stress_value
        );

        const Index row = static_cast<Index>(offset) + point;

        if (strain) {
            for (Index component = 0; component < strain->components; ++component) {
                (*strain)(row, component) = component < 6
                    ? strain_value(component)
                    : Precision(0);
            }
        }

        if (stress) {
            for (Index component = 0; component < stress->components; ++component) {
                (*stress)(row, component) = component < 6
                    ? stress_value(component)
                    : Precision(0);
            }
        }
    }
}

/**
 * Computes the eight generalized shell resultants at all integration points.
 *
 * @param stress_state Output integration-point resultant field.
 * @param displacement Global nodal displacement field.
 * @param offset First output row belonging to this element.
 * @param use_green_lagrange_nl Select nonlinear or linearized evaluation.
 */
template<Index N>
void FRTShell<N>::compute_stress_state(Field&       stress_state,
                                       const Field& displacement,
                                       int          offset,
                                       bool         use_green_lagrange_nl) {
    logging::error(stress_state.components >= num_strains,
                   "FRTShell: stress state requires at least eight components "
                   "[N11,N22,N12,M11,M22,M12,Q13,Q23]");

    if (use_green_lagrange_nl) {
        const CurrentState state = current_state_from_displacement(displacement);
        const EvaluationData data = init_evaluation(
            state,
            true,
            false,
            false,
            true
        );

        for (Index ip = 0; ip < static_cast<Index>(data.ip_resultants.size()); ++ip) {
            const Index row = static_cast<Index>(offset) + ip;

            for (Index component = 0; component < stress_state.components; ++component) {
                stress_state(row, component) = component < num_strains
                    ? data.ip_resultants[static_cast<std::size_t>(ip)](component)
                    : Precision(0);
            }
        }

        return;
    }

    const CurrentState state = reference_state();
    const EvaluationData data = init_evaluation(
        state,
        true,
        true,
        false,
        false
    );
    const Vec6N q = element_displacement_vector(displacement);

    for (Index ip = 0; ip < static_cast<Index>(data.ip_B.size()); ++ip) {
        const std::size_t id     = static_cast<std::size_t>(ip);
        const Vec8        strain = data.ip_B[id] * q;
        const Vec8        values = data.ip_tangent[id] * strain;
        const Index       row    = static_cast<Index>(offset) + ip;

        for (Index component = 0; component < stress_state.components; ++component) {
            stress_state(row, component) = component < num_strains
                ? values(component)
                : Precision(0);
        }
    }
}

/**
 * Averages generalized shell resultants from the element nodes into nodal
 * result fields.
 *
 * @param resultants Global nodal generalized-resultant accumulator.
 * @param contribution_count Global nodal contribution counter.
 * @param displacement Global nodal displacement field.
 * @return Always `true` after shell resultants were accumulated.
 */
template<Index N>
bool FRTShell<N>::compute_shell_section_forces(Field&       resultants,
                                               Field&       contribution_count,
                                               const Field& displacement) {
    logging::error(resultants.components >= num_strains,
                   "FRTShell: shell section forces require eight components "
                   "[N11,N22,N12,M11,M22,M12,Q13,Q23]");

    const RowMatrix rst = this->stress_strain_nodal_rst();
    const CurrentState state = current_state_from_displacement(displacement);
    const EvaluationData data = init_evaluation(
        state,
        true,
        false,
        false,
        false
    );
    const Vec6N q = element_displacement_vector(displacement);

    for (Index node = 0; node < num_nodes; ++node) {
        const Vec8 values = generalized_resultant_at(
            data,
            q,
            rst(node, 0),
            rst(node, 1),
            true
        );
        const Index node_id = static_cast<Index>(this->node_ids[node]);

        for (Index component = 0; component < num_strains; ++component) {
            resultants(node_id, component) += values(component);
        }

        contribution_count(node_id, 0) += Precision(1);
    }

    return true;
}

} // namespace fem::model
