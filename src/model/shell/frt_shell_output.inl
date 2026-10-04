/**
 * @file frt_shell_output.inl
 * @brief Implements generalized and physical shell result recovery.
 *
 * The routines evaluate MITC generalized strains and shell resultants at
 * arbitrary natural coordinates. Physical through-thickness strains are
 * reconstructed from generalized shell strains. Material stresses are evaluated
 * by the shell section and returned in the configured stress-output convention.
 *
 * Output recovery is state-neutral: committed history may be inspected, but no
 * result query receives a writable persistent trial-state target.
 *
 * @see FRTShell
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#include "frt_shell.h"

#include "../../core/logging.h"
#include "../../material/isotropic_j2_elasticity.h"
#include "../../material/strain/volume_strain_green_lagrange.h"
#include "../../material/stress/volume_stress_cauchy.h"
#include "../../math/extrapolate.h"
#include "../../math/vec_util.h"

#include <algorithm>
#include <cmath>

namespace fem::model {

using math::normalized;

template<Index N>
Precision FRTShell<N>::thermal_free_strain_at(const Field* thermal_free_strain,
                                              Precision r,
                                              Precision s) const {
    if (!thermal_free_strain) {
        return Precision(0);
    }

    logging::error(thermal_free_strain->domain == FieldDomain::ELEMENT_NODAL
                   && thermal_free_strain->components == 1,
                   "FRTShell: thermal free strain must be scalar ELEMENT_NODAL data");

    const VecN shape = shape_function(r, s);
    Precision value = Precision(0);
    for (Index node = 0; node < num_nodes; ++node) {
        value += shape(node)
            * (*thermal_free_strain)(
                static_cast<Index>(this->elem_nodal_offset) + node, 0);
    }
    return value;
}

/**
 * Evaluates generalized shell strain at an arbitrary natural point by
 * linearizing the exact finite shell kinematics about the state stored in data.
 *
 *     epsilon ~= epsilon0 + B0 Delta q.
 *
 * Passing a zero increment therefore recovers the exact base strain.
 */
template<Index N>
typename FRTShell<N>::Vec8 FRTShell<N>::generalized_strain_at(
    const EvaluationData& data,
    const Vec6N&          displacement_increment,
    Precision             r,
    Precision             s
) const {
    logging::error(data.with_B,
        "FRTShell: affine generalized strain recovery requires base-state B matrices");

    const ReferencePoint* cached    = cached_reference_point(r, s);
    const ReferencePoint  temporary = cached ? ReferencePoint{} : make_reference_point(r, s, Precision(0));
    const ReferencePoint& point     = cached ? *cached : temporary;

    Vec8     strain;
    Mat8x6N B;

    compute_natural_strain(data, point, strain, &B);
    apply_mitc_natural(data, point, strain, &B);
    transform_strain_to_local(point, strain, &B);

    return strain + B * displacement_increment;
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

        x_a += point.shape_ab.col(0)(node) * x_i;
        x_b += point.shape_ab.col(1)(node) * x_i;
        d   += point.shape(node) * d_i;
        d_a += point.shape_ab.col(0)(node) * d_i;
        d_b += point.shape_ab.col(1)(node) * d_i;
    }

    Mat3 reference_covariant;
    reference_covariant.col(0) = point.X_ab.col(0) + z * point.D_ab.col(0);
    reference_covariant.col(1) = point.X_ab.col(1) + z * point.D_ab.col(1);
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
 * Returns the first variation of the three-dimensional shell deformation
 * gradient about the base state stored in data.
 *
 * Translational increments vary the midsurface tangents directly. Rotational
 * increments use the exact first SO(3) derivatives at q0:
 *
 *     Delta d_i = sum_a dR_i/dtheta_a d0_i Delta theta_a.
 */
template<Index N>
Mat3 FRTShell<N>::deformation_gradient_increment_at(
    const EvaluationData& data,
    const Vec6N&          displacement_increment,
    Precision             r,
    Precision             s,
    Precision             z
) const {
    logging::error(data.rotations != nullptr,
        "FRTShell: deformation-gradient linearization requires base-state rotation derivatives");

    const ReferencePoint point = make_reference_point(r, s, Precision(0));
    const auto& rotations = *data.rotations;
    const auto& d0        = reference_data().d0;

    Vec3 delta_x_a = Vec3::Zero();
    Vec3 delta_x_b = Vec3::Zero();
    Vec3 delta_d   = Vec3::Zero();
    Vec3 delta_d_a = Vec3::Zero();
    Vec3 delta_d_b = Vec3::Zero();

    for (Index node = 0; node < num_nodes; ++node) {
        const Vec3 delta_x     = displacement_increment.template segment<3>(dofs_per_node * node);
        const Vec3 delta_theta = displacement_increment.template segment<3>(dofs_per_node * node + 3);
        const Vec3 d0_i        = d0.row(node).transpose();

        Vec3 delta_director = Vec3::Zero();
        for (Index component = 0; component < 3; ++component) {
            delta_director.noalias() += delta_theta(component) * (rotations[node].d1[component] * d0_i);
        }

        delta_x_a += point.shape_ab(node, 0) * delta_x;
        delta_x_b += point.shape_ab(node, 1) * delta_x;

        delta_d   += point.shape(node)       * delta_director;
        delta_d_a += point.shape_ab(node, 0) * delta_director;
        delta_d_b += point.shape_ab(node, 1) * delta_director;
    }

    Mat3 reference_covariant;
    reference_covariant.col(0) = point.X_ab.col(0) + z * point.D_ab.col(0);
    reference_covariant.col(1) = point.X_ab.col(1) + z * point.D_ab.col(1);
    reference_covariant.col(2) = point.D;

    const Precision reference_det = reference_covariant.determinant();
    logging::error(std::abs(reference_det) > Precision(1e-14),
        "FRTShell: singular reference shell basis during stress linearization in element ", this->elem_id);

    Mat3 delta_covariant;
    delta_covariant.col(0) = delta_x_a + z * delta_d_a;
    delta_covariant.col(1) = delta_x_b + z * delta_d_b;
    delta_covariant.col(2) = delta_d;

    return delta_covariant * reference_covariant.inverse();
}

/**
 * Reconstructs physical through-thickness strain and stress tensors from an
 * exact base state followed by a first-order perturbation.
 *
 * Generalized strain is continued as
 *
 *     epsilon ~= epsilon0 + B0 Delta q.
 *
 * The section evaluates PK2 stress and its tangent at epsilon0. Cauchy stress is
 * then linearized consistently with respect to both Delta S and Delta F.
 */
template<Index N>
void FRTShell<N>::physical_stress_strain_at(
    const EvaluationData& data,
    const Vec6N&          displacement_increment,
    Precision             r,
    Precision             s,
    Precision             zeta,
    const Field*          thermal_free_strain,
    Vec6&                 strain_out,
    Vec6&                 stress_out
) const {
    const Vec8 generalized_base   = generalized_strain_at(data, Vec6N::Zero(), r, s);
    const Vec8 generalized_strain = generalized_strain_at(data, displacement_increment, r, s);
    Vec8 generalized_increment    = generalized_strain - generalized_base;

    if (thermal_free_strain) {
        const ReferencePoint* cached    = cached_reference_point(r, s);
        const ReferencePoint  temporary = cached ? ReferencePoint{} : make_reference_point(r, s, Precision(0));
        const ReferencePoint& point     = cached ? *cached : temporary;
        generalized_increment -= thermal_generalized_strain(point, thermal_free_strain_at(thermal_free_strain, r, s));
    }

    const Precision h = this->get_section()->thickness_;
    const Precision z = Precision(0.5) * h * zeta;

    const Index membrane_strain_start = static_cast<Index>(ShellGeneralizedStrain::Component::EpsilonXX);
    const Index curvature_start       = static_cast<Index>(ShellGeneralizedStrain::Component::KappaXX);
    const Index shear_strain_start    = static_cast<Index>(ShellGeneralizedStrain::Component::GammaXZ);

    const Vec3 plane_strain = generalized_strain.template segment<3>(membrane_strain_start)
                            + z * generalized_strain.template segment<3>(curvature_start);
    const Vec2 shear_strain = generalized_strain.template segment<2>(shear_strain_start);

    VolumeStrainGreenLagrange strain_local;
    strain_local[VolumeStrain::Component::XX]      = plane_strain(0);
    strain_local[VolumeStrain::Component::YY]      = plane_strain(1);
    strain_local[VolumeStrain::Component::GammaYZ] = shear_strain(1);
    strain_local[VolumeStrain::Component::GammaXZ] = shear_strain(0);
    strain_local[VolumeStrain::Component::GammaXY] = plane_strain(2);

    const Mat3 reference_basis       = reference_basis_global(r, s);
    const Mat3 green_lagrange_global = reference_basis * strain_local.tensor() * reference_basis.transpose();
    strain_out = VolumeStrainGreenLagrange(green_lagrange_global).voigt();

    const Mat3 deformation_gradient_base = deformation_gradient_at(data.state, r, s, z);
    const Mat3 deformation_gradient_increment =
        deformation_gradient_increment_at(data, displacement_increment, r, s, z);

    const auto& points = reference_data().ip_points;
    Index state_ip = 0;
    Precision state_distance = (r - points[0].r) * (r - points[0].r) + (s - points[0].s) * (s - points[0].s);

    for (Index ip = 1; ip < static_cast<Index>(points.size()); ++ip) {
        const ReferencePoint& point = points[static_cast<std::size_t>(ip)];
        const Precision distance = (r - point.r) * (r - point.r) + (s - point.s) * (s - point.s);
        if (distance < state_distance) {
            state_ip       = ip;
            state_distance = distance;
        }
    }

    const Index      state_row = this->mp_index(state_ip, 0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);

    const VolumeStressCauchy cauchy_stress = this->get_section()->recover_stress(
        reference_position(r, s),
        reference_basis,
        ShellGeneralizedStrain(generalized_base),
        ShellGeneralizedStrain(generalized_increment),
        old_state,
        this->_model_data->material_state_old->components,
        z,
        deformation_gradient_base,
        deformation_gradient_increment
    );

    stress_out = topology_stiffness_scale() * cauchy_stress.voigt();
}

/**
 * Computes physical stress and strain at displacement q by linearizing about
 * the exact finite-rotation state q0 supplied through linearization.
 *
 * A null linearization denotes q0 = 0. Passing displacement itself recovers the
 * exact nonlinear state because Delta q then vanishes.
 */
template<Index N>
void FRTShell<N>::compute_stress_strain(
    Field*           strain,
    Field*           stress,
    const Field&     displacement,
    const RowMatrix& rst,
    int              offset,
    const Field*     linearization,
    const Field*     thermal_free_strain
) {
    logging::error(strain != nullptr || stress != nullptr,
        "FRTShell: compute_stress_strain requires at least one output field");
    logging::error(rst.cols() >= 3,
        "FRTShell: stress/strain coordinates require r, s and t columns");
    logging::error(!thermal_free_strain || linearization == nullptr,
        "FRTShell: finite-state thermal free strain recovery is not implemented");
    logging::error(!thermal_free_strain || (thermal_free_strain->domain == FieldDomain::ELEMENT_NODAL && thermal_free_strain->components == 1),
        "FRTShell: thermal free strain must be scalar ELEMENT_NODAL data");

    const Vec6N q0 = linearization ? element_displacement_vector(*linearization) : Vec6N::Zero();
    const Vec6N q  = element_displacement_vector(displacement);
    const Vec6N delta = q - q0;

    const CurrentState state = linearization ? current_state_from_displacement(*linearization) : reference_state();

    // Request B0 and the first SO(3) derivatives at the exact base state. The
    // resultant request keeps section evaluation on the same Green-Lagrange
    // constitutive path used by mechanical evaluation.
    const EvaluationData data = init_evaluation(state, true, true, false, true, false);

    for (Eigen::Index point = 0; point < rst.rows(); ++point) {
        Vec6 strain_value;
        Vec6 stress_value;

        physical_stress_strain_at(data, delta, rst(point, 0), rst(point, 1), rst(point, 2),
                                  thermal_free_strain, strain_value, stress_value);

        const Index row = static_cast<Index>(offset) + point;

        if (strain) {
            for (Index component = 0; component < strain->components; ++component) {
                (*strain)(row, component) = component < 6 ? strain_value(component) : Precision(0);
            }
        }

        if (stress) {
            for (Index component = 0; component < stress->components; ++component) {
                (*stress)(row, component) = component < 6 ? stress_value(component) : Precision(0);
            }
        }
    }
}

/**
 * Recovers accumulated equivalent plastic strain from committed shell material
 * history.
 *
 * Every in-plane shell integration point owns a contiguous block of
 * through-thickness material states. The largest J2 PEEQ value in that block is
 * retained as the scalar shell value for the in-plane point. Those maxima are
 * then extrapolated in natural midsurface coordinates to the shell nodes.
 *
 * Linear triangular recovery is used for S3 and S6 because their common cubic
 * triangle rule does not provide enough independent samples for a complete
 * quadratic field. S4 uses bilinear recovery and S8 uses the quadratic
 * serendipity basis supported by its 3 x 3 integration rule.
 *
 * A J2 material whose nonlinear state storage has not been initialized yet
 * contributes zero PEEQ. Non-J2 shell sections do not participate in the
 * model-wide PEEQ average.
 *
 * @param peeq Scalar ELEMENT_NODAL output field.
 * @param offset First element-nodal row belonging to this shell.
 * @return True when the shell material uses J2 plasticity.
 */
template<Index N>
bool FRTShell<N>::compute_peeq(Field& peeq, int offset) {
    logging::error(peeq.domain == FieldDomain::ELEMENT_NODAL && peeq.components == 1,
        "FRTShell: PEEQ recovery requires scalar ELEMENT_NODAL output");

    auto mat = this->get_material();
    if (!mat || !mat->has_elasticity()) return false;

    const auto* j2 = mat->elasticity()->template as<material::IsotropicJ2Elasticity>();
    if (!j2) return false;

    const RowMatrix ip_rst    = this->stress_strain_ip_rst();
    const RowMatrix nodal_rst = this->stress_strain_nodal_rst();
    RowMatrix ip_peeq         = RowMatrix::Zero(ip_rst.rows(), 1);

    // A virgin J2 shell has zero accumulated plastic strain before nonlinear
    // state storage has been expanded to the constitutive state width.
    const auto& state = this->_model_data->material_state_old;
    if (state && state->components >= j2->state_size()) {
        const Index material_points = this->num_mp_per_ip();

        // Reduce the through-thickness history to its maximum at each in-plane IP.
        for (Index ip = 0; ip < static_cast<Index>(ip_rst.rows()); ++ip) {
            Precision maximum = Precision(0);

            for (Index mp = 0; mp < material_points; ++mp) {
                const Index state_row = this->mp_index(ip, mp);
                const Precision* old_state = &(*state)(state_row, 0);
                maximum = std::max(maximum, j2->equivalent_plastic_strain(old_state));
            }

            ip_peeq(ip, 0) = maximum;
        }
    }

    // Reconstruct the reduced scalar field at the natural shell nodes. The
    // operator depends only on topology and quadrature, so build it once per
    // concrete FRT shell topology instead of refactorizing it for every element.
    using math::ExtrapolationBasis;
    RowMatrix nodal_peeq;

    if constexpr (N == 3 || N == 6) {
        static const RowMatrix E = math::extrapolate(
            ip_rst,
            nodal_rst,
            {ExtrapolationBasis::F1, ExtrapolationBasis::FR, ExtrapolationBasis::FS}
        );
        nodal_peeq = E * ip_peeq;
    } else if constexpr (N == 4) {
        static const RowMatrix E = math::extrapolate(
            ip_rst,
            nodal_rst,
            {
                ExtrapolationBasis::F1,
                ExtrapolationBasis::FR,
                ExtrapolationBasis::FS,
                ExtrapolationBasis::FRS
            }
        );
        nodal_peeq = E * ip_peeq;
    } else {
        static const RowMatrix E = math::extrapolate(
            ip_rst,
            nodal_rst,
            {
                ExtrapolationBasis::F1,
                ExtrapolationBasis::FR,
                ExtrapolationBasis::FS,
                ExtrapolationBasis::FRR,
                ExtrapolationBasis::FSS,
                ExtrapolationBasis::FRS,
                ExtrapolationBasis::FRRS,
                ExtrapolationBasis::FSSR
            }
        );
        nodal_peeq = E * ip_peeq;
    }
    for (Index node = 0; node < N; ++node) {
        peeq(static_cast<Index>(offset) + node, 0) = nodal_peeq(node, 0);
    }

    return true;
}

/**
 * Averages generalized shell resultants from the element nodes into nodal
 * result fields without advancing material trial history.
 *
 * Each natural nodal output coordinate reuses the closest in-plane integration-
 * point committed state block. This recovery is used by the linear-static load
 * case and therefore evaluates linearized generalized strains in the reference
 * shell configuration. The section evaluates resultants in its configured
 * output basis before they are accumulated for subsequent component-wise
 * averaging.
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
    return compute_shell_section_forces(
        resultants, contribution_count, displacement, nullptr);
}

template<Index N>
bool FRTShell<N>::compute_shell_section_forces(Field&       resultants,
                                               Field&       contribution_count,
                                               const Field& displacement,
                                               const Field* thermal_free_strain) {
    logging::error(resultants.components >= num_strains,
        "FRTShell: shell section forces require eight components [N11,N22,N12,M11,M22,M12,Q13,Q23]");
    logging::error(!thermal_free_strain || (thermal_free_strain->domain == FieldDomain::ELEMENT_NODAL && thermal_free_strain->components == 1),
        "FRTShell: thermal free strain must be scalar ELEMENT_NODAL data");

    const RowMatrix      rst     = this->stress_strain_nodal_rst();
    const CurrentState   state   = reference_state();
    const EvaluationData data    = init_evaluation(state, true, true, false, false);
    const Vec6N          q       = element_displacement_vector(displacement);
    ShellSection*        section = shell_section();
    const Precision      scale   = topology_stiffness_scale();
    const auto&          points  = reference_data().ip_points;

    // Recover one generalized resultant vector at every natural element node.
    //
    // The reference state q0 = 0 is evaluated exactly and the requested linear
    // displacement state is continued affinely:
    //
    //     n(q) ~= n0 + H0 Delta epsilon.
    for (Index node = 0; node < num_nodes; ++node) {
        const Precision r = rst(node, 0);
        const Precision s = rst(node, 1);

        const Vec8 strain_base = generalized_strain_at(data, Vec6N::Zero(), r, s);
        Vec8 strain_increment  = generalized_strain_at(data, q, r, s) - strain_base;

        if (thermal_free_strain) {
            const ReferencePoint* cached    = cached_reference_point(r, s);
            const ReferencePoint  temporary = cached ? ReferencePoint{} : make_reference_point(r, s, Precision(0));
            const ReferencePoint& point     = cached ? *cached : temporary;
            strain_increment -= thermal_generalized_strain(point, thermal_free_strain_at(thermal_free_strain, r, s));
        }

        // Natural nodal output points own no independent constitutive history.
        // Reuse the committed block of the closest in-plane integration point.
        Index state_ip = 0;
        Precision state_distance =
            (r - points[0].r) * (r - points[0].r)
            + (s - points[0].s) * (s - points[0].s);

        for (Index ip = 1; ip < static_cast<Index>(points.size()); ++ip) {
            const ReferencePoint& point = points[static_cast<std::size_t>(ip)];
            const Precision distance =
                (r - point.r) * (r - point.r)
                + (s - point.s) * (s - point.s);

            if (distance < state_distance) {
                state_ip       = ip;
                state_distance = distance;
            }
        }

        const Index      state_row = this->mp_index(state_ip, 0);
        const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);

        ShellStressResultants resultants_base;
        Mat8                  tangent_base;

        const Vec3 position    = reference_position(r, s);
        const Mat3 shell_basis = reference_basis_global(r, s);

        section->evaluate(position, shell_basis, ShellGeneralizedStrain(strain_base), old_state, nullptr, this->_model_data->material_state_old->components, resultants_base, tangent_base);

        ShellStressResultants resultants_shell(
            resultants_base.values() + tangent_base * strain_increment
        );

        // Resultants from neighboring shells must be expressed in one
        // deterministic tangential basis before component-wise nodal averaging.
        // This is an element-output convention, not a constitutive section task.
        const Precision projection_tolerance = Precision(1e-6);
        const Vec3      normal               = shell_basis.col(2).normalized();

        Vec3 source_axis;
        if (section->orientation_) {
            const Vec3 point_local      = section->orientation_->to_local(position);
            const Mat3 orientation_axes = section->orientation_->get_axes(point_local);
            source_axis = orientation_axes.col(section->csys_axis_);
        } else {
            source_axis = Vec3::UnitX();
        }

        Vec3 resultant_e1 = source_axis - normal * source_axis.dot(normal);

        if (!section->orientation_ && resultant_e1.norm() <= projection_tolerance) {
            source_axis  = Vec3::UnitY();
            resultant_e1 = source_axis - normal * source_axis.dot(normal);
        }

        logging::error(resultant_e1.norm() > projection_tolerance,
            "FRTShell: selected output axis cannot define a tangential resultant basis for element ", this->elem_id);

        resultant_e1.normalize();
        const Vec3 resultant_e2 = normal.cross(resultant_e1).normalized();

        Mat3 result_basis;
        result_basis.col(0) = resultant_e1;
        result_basis.col(1) = resultant_e2;
        result_basis.col(2) = normal;

        const Mat2 result_axes_in_shell =
            shell_basis.template block<3, 2>(0, 0).transpose()
            * result_basis.template block<3, 2>(0, 0);

        resultants_shell = resultants_shell.transformed(result_axes_in_shell);
        const Vec8 values = scale * resultants_shell.values();

        const Index node_id = static_cast<Index>(this->node_ids[node]);
        for (Index component = 0; component < num_strains; ++component) {
            resultants(node_id, component) += values(component);
        }

        contribution_count(node_id, 0) += Precision(1);
    }

    return true;
}

} // namespace fem::model
