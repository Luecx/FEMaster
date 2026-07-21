/**
 * @file frt_shell_kinematics.inl
 * @brief Implements finite-rotation shell kinematics and compact derivatives.
 *
 * The implementation evaluates the current nodal directors through the SO(3)
 * exponential map and constructs compatible Total-Lagrangian membrane,
 * bending and transverse-shear strains in natural coordinates. First strain
 * derivatives are assembled directly from analytical translational identities
 * and compact local rotation derivatives.
 *
 * No element-sized position derivative matrices and no dense generalized
 * strain Hessians are stored. Second SO(3) derivatives are retained only in the
 * active thread-local workspace when the geometric tangent is requested.
 *
 * @see FRTShell
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#include "frt_shell.h"

#include "../../core/logging.h"
#include "../../math/so3.h"

#include <cmath>

namespace fem::model {

namespace {


/**
 * Transforms one three-component engineering tensor block and its optional B
 * rows by a constant pointwise reference transformation.
 *
 * The helper is used for both the membrane and curvature blocks. The source is
 * copied before assignment so in-place transformations remain valid.
 *
 * @tparam N Number of shell nodes.
 * @param transformation Engineering tensor transformation matrix.
 * @param start First component of the three-row block.
 * @param strain Generalized strain vector to transform in place.
 * @param B Optional generalized strain-displacement matrix to transform.
 */
template<Index N>
void transform_in_plane_rows(
    const StaticMatrix<3, 3>&      transformation,
    Index                          start,
    typename FRTShell<N>::Vec8&    strain,
    typename FRTShell<N>::Mat8x6N* B
) {
    using Shell = FRTShell<N>;

    // Copy the current values because source and target occupy the same block
    const StaticVector<3> values = strain.template segment<3>(start);
    strain.template segment<3>(start) = transformation * values;

    if (!B) {
        return;
    }

    // Apply the identical linear map to the corresponding B rows
    const StaticMatrix<3, Shell::num_dofs> rows =
        B->template block<3, Shell::num_dofs>(start, 0);
    B->template block<3, Shell::num_dofs>(start, 0) = transformation * rows;
}

} // namespace

/**
 * Collects the current global nodal positions and total rotation vectors.
 *
 * The model POSITION field stores six values per shell node in the order
 * `[x, y, z, rx, ry, rz]` for the current nonlinear configuration.
 *
 * @return Matrix containing one current six-component nodal state per row.
 */
template<Index N>
typename FRTShell<N>::MatN6 FRTShell<N>::node_coords_current_6() const {
    logging::error(this->_model_data            != nullptr,
                   "FRTShell: no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions != nullptr,
                   "FRTShell: POSITION field is not set");

    const auto& positions = *this->_model_data->positions;
    MatN6      values;

    // Gather the element rows once in local nodal ordering
    for (Index node = 0; node < num_nodes; ++node) {
        values.row(node) = positions.row_vec6(
            static_cast<Index>(this->node_ids[node])
        ).transpose();
    }

    return values;
}

/**
 * Constructs the current nodal shell state from the model POSITION field.
 *
 * The total nodal axis-angle vector is converted to one SO(3) rotation matrix
 * and applied to the corresponding reference director.
 *
 * @return Current nodal positions, directors and total rotations.
 */
template<Index N>
typename FRTShell<N>::CurrentState FRTShell<N>::current_state() const {
    const auto& ref       = reference_data();
    const MatN6 positions = node_coords_current_6();

    CurrentState state;

    // Construct the physical nodal position and director at every shell node
    for (Index node = 0; node < num_nodes; ++node) {
        const Vec3 x     = positions.template block<1, 3>(node, 0).transpose();
        const Vec3 theta = positions.template block<1, 3>(node, 3).transpose();
        const Vec3 d0    = ref.d0.row(node).transpose();
        const Mat3 R     = math::so3::rotation_matrix(theta);

        state.x.row(node)     = x.transpose();
        state.d.row(node)     = (R * d0).transpose();
        state.theta.row(node) = theta.transpose();
    }

    return state;
}

/**
 * Constructs the undeformed nodal state used by linear shell evaluation.
 *
 * @return Reference positions and directors with zero total rotations.
 */
template<Index N>
typename FRTShell<N>::CurrentState FRTShell<N>::reference_state() const {
    const auto& ref = reference_data();

    CurrentState state;
    state.x = ref.X;
    state.d = ref.d0;
    state.theta.setZero();
    return state;
}

/**
 * Constructs a trial shell state from the reference geometry and a supplied
 * nodal displacement field.
 *
 * Translational components are added to the reference midsurface coordinates.
 * Rotational components are interpreted as total global axis-angle vectors and
 * rotate the nodal reference directors through the SO(3) exponential map.
 *
 * @param displacement Nodal six-component trial displacement field.
 * @return Trial nodal positions, directors and total rotations.
 */
template<Index N>
typename FRTShell<N>::CurrentState FRTShell<N>::current_state_from_displacement(
    const Field& displacement
) const {
    logging::error(displacement.domain == FieldDomain::NODE,
                   "FRTShell: displacement field must use NODE domain");
    logging::error(displacement.components >= dofs_per_node,
                   "FRTShell: displacement field requires six components");

    const auto& ref = reference_data();
    CurrentState state;

    // Build the complete trial state in local element-node ordering
    for (Index node = 0; node < num_nodes; ++node) {
        const Index node_id = static_cast<Index>(this->node_ids[node]);
        const Vec6  q       = displacement.row_vec6(node_id);
        const Vec3  x       = ref.X.row(node).transpose() + q.template head<3>();
        const Vec3  theta   = q.template tail<3>();
        const Vec3  d0      = ref.d0.row(node).transpose();
        const Mat3  R       = math::so3::rotation_matrix(theta);

        state.x.row(node)     = x.transpose();
        state.d.row(node)     = (R * d0).transpose();
        state.theta.row(node) = theta.transpose();
    }

    return state;
}

/**
 * Gathers the element displacement vector in nodal six-degree-of-freedom
 * ordering.
 *
 * @param displacement Global nodal displacement field.
 * @return Fixed-size element displacement vector.
 */
template<Index N>
typename FRTShell<N>::Vec6N FRTShell<N>::element_displacement_vector(
    const Field& displacement
) const {
    logging::error(displacement.domain == FieldDomain::NODE,
                   "FRTShell: displacement field must use NODE domain");
    logging::error(displacement.components >= dofs_per_node,
                   "FRTShell: displacement field requires six components");

    Vec6N q;

    // Gather the six nodal components without intermediate dynamic storage
    for (Index node = 0; node < num_nodes; ++node) {
        const Index node_id = static_cast<Index>(this->node_ids[node]);
        q.template segment<6>(dofs_per_node * node) = displacement.row_vec6(node_id);
    }

    return q;
}

/**
 * Evaluates the compatible generalized shell strain in natural coordinates.
 *
 * The generalized natural strain vector is ordered as
 *
 *     strain_nat =
 *     [ epsilon_rr,
 *       epsilon_ss,
 *       gamma_rs,
 *       kappa_rr,
 *       kappa_ss,
 *       kappa_rs,
 *       gamma_r3,
 *       gamma_s3 ].
 *
 * Membrane strains are obtained from changes of the current midsurface metric.
 * Bending strains are obtained from products between current midsurface
 * tangents and current director derivatives. Transverse-shear strains compare
 * the current midsurface tangents with the interpolated current director.
 *
 * The current contributions are evaluated first. The corresponding reference
 * metric, curvature and shear quantities are subtracted afterwards, such that
 * the complete generalized strain vanishes in the undeformed configuration.
 *
 * When `B_nat` is supplied, it contains the first derivative
 *
 *     B_nat = d strain_nat / d q_e,
 *
 * where the element degrees of freedom are arranged node-wise as
 *
 *     q_i = [u_x,i, u_y,i, u_z,i, theta_x,i, theta_y,i, theta_z,i].
 *
 * Each helper inserts its analytical derivative directly into the row belonging
 * to the generalized strain component being evaluated. Translational
 * derivatives are assembled through compact matrix products and strided Eigen
 * indexing. Director derivatives use the cached nodal SO(3) derivatives.
 *
 * @param data Active element evaluation data.
 * @param point Pointwise curved-reference geometry.
 * @param strain_nat Output compatible generalized strain in natural components.
 * @param B_nat Optional output first strain derivative.
 */
template<Index N>
void FRTShell<N>::compute_natural_strain(
    const EvaluationData& data,
    const ReferencePoint& point,
    Vec8&                 strain_nat,
    Mat8x6N*              B_nat
) const {
    using Component = ShellGeneralizedStrain::Component;

    constexpr Index epsilon_rr = static_cast<Index>(Component::EpsilonXX);
    constexpr Index epsilon_ss = static_cast<Index>(Component::EpsilonYY);
    constexpr Index gamma_rs   = static_cast<Index>(Component::GammaXY);
    constexpr Index kappa_rr   = static_cast<Index>(Component::KappaXX);
    constexpr Index kappa_ss   = static_cast<Index>(Component::KappaYY);
    constexpr Index kappa_rs   = static_cast<Index>(Component::KappaXY);
    constexpr Index gamma_r3   = static_cast<Index>(Component::GammaXZ);
    constexpr Index gamma_s3   = static_cast<Index>(Component::GammaYZ);

    strain_nat.setZero();

    if (B_nat) {
        B_nat->setZero();
    }

    // Store the pointwise shape-function values and their natural derivatives
    // once because they are reused by the membrane, bending and transverse-
    // shear strain terms:
    //
    //     shape_r(i) = N_i,r,
    //     shape_s(i) = N_i,s,
    //     shape(i)   = N_i.
    const VecN shape_r = point.shape_rs.col(0);
    const VecN shape_s = point.shape_rs.col(1);
    const VecN shape   = point.shape;
    const MatN3& d0    = reference_data().d0;

    // -------------------------------------------------------------------------
    // Membrane strains
    // -------------------------------------------------------------------------
    //
    // The current midsurface position is interpolated as
    //
    //     x(r,s) = sum_i N_i(r,s) x_i.
    //
    // Its natural tangent vectors are
    //
    //     x_,r = sum_i N_i,r x_i,
    //     x_,s = sum_i N_i,s x_i.
    //
    // The current covariant metric components are
    //
    //     g_rr = dot(x_,r, x_,r),
    //     g_ss = dot(x_,s, x_,s),
    //     g_rs = dot(x_,r, x_,s).
    //
    // The Green-Lagrange membrane components use
    //
    //     epsilon_rr,current = 0.5 g_rr,
    //     epsilon_ss,current = 0.5 g_ss,
    //     gamma_rs,current   =       g_rs.
    //
    // `gamma_rs` is stored as an engineering shear component and therefore does
    // not contain the factor 0.5.
    //
    // evaluate_xx_product() returns the corresponding scalar product and,
    // when B_nat is supplied, inserts its derivative directly into the
    // translational columns of the selected B_nat row. The rotational columns
    // of these three rows remain zero because the membrane metric depends only
    // on the current nodal positions.

    strain_nat(epsilon_rr) += evaluate_xx_product<N>(
        data.state.x, shape_r, shape_r, Precision(0.5), epsilon_rr, B_nat
    );
    strain_nat(epsilon_ss) += evaluate_xx_product<N>(
        data.state.x, shape_s, shape_s, Precision(0.5), epsilon_ss, B_nat
    );
    strain_nat(gamma_rs) += evaluate_xx_product<N>(
        data.state.x, shape_r, shape_s, Precision(1), gamma_rs, B_nat
    );

    // -------------------------------------------------------------------------
    // Bending strains
    // -------------------------------------------------------------------------
    //
    // The current nodal directors are
    //
    //     d_i = R_i d0_i,
    //
    // and the current interpolated director field is
    //
    //     d(r,s) = sum_i N_i d_i.
    //
    // Its natural derivatives are
    //
    //     d_,r = sum_i N_i,r d_i,
    //     d_,s = sum_i N_i,s d_i.
    //
    // The current curvature terms are formed from products between the current
    // midsurface tangents and current director derivatives:
    //
    //     kappa_rr,current = dot(x_,r, d_,r),
    //     kappa_ss,current = dot(x_,s, d_,s),
    //
    //     kappa_rs,current =
    //         dot(x_,r, d_,s) + dot(x_,s, d_,r).
    //
    // The mixed curvature component is stored in engineering convention and is
    // therefore assembled from both off-diagonal contributions.
    //
    // evaluate_xd_product() inserts both derivative
    // parts directly into the selected B_nat row:
    //
    //   - translation derivatives caused by x_,r or x_,s,
    //   - rotation derivatives caused by d_,r or d_,s.

    strain_nat(kappa_rr) += evaluate_xd_product<N>(
        data, d0, shape_r, shape_r, Precision(1), kappa_rr, B_nat
    );
    strain_nat(kappa_ss) += evaluate_xd_product<N>(
        data, d0, shape_s, shape_s, Precision(1), kappa_ss, B_nat
    );
    strain_nat(kappa_rs) += evaluate_xd_product<N>(
        data, d0, shape_r, shape_s, Precision(1), kappa_rs, B_nat
    );
    strain_nat(kappa_rs) += evaluate_xd_product<N>(
        data, d0, shape_s, shape_r, Precision(1), kappa_rs, B_nat
    );

    // -------------------------------------------------------------------------
    // Transverse-shear strains
    // -------------------------------------------------------------------------
    //
    // The interpolated current director at the evaluation point is
    //
    //     d = sum_i N_i d_i.
    //
    // The two covariant transverse-shear measures are
    //
    //     gamma_r3,current = dot(x_,r, d),
    //     gamma_s3,current = dot(x_,s, d).
    //
    // The position coefficients are therefore the natural shape derivatives,
    // while the director coefficients are the shape-function values
    // themselves.

    strain_nat(gamma_r3) += evaluate_xd_product<N>(
        data, d0, shape_r, shape, Precision(1), gamma_r3, B_nat
    );

    strain_nat(gamma_s3) += evaluate_xd_product<N>(
        data, d0, shape_s, shape, Precision(1), gamma_s3, B_nat
    );

    // -------------------------------------------------------------------------
    // Reference-state subtraction
    // -------------------------------------------------------------------------
    //
    // The terms evaluated above contain the complete current metric, curvature
    // and transverse-shear quantities. Subtracting their pointwise reference
    // values produces the corresponding changes:
    //
    //     epsilon_rr = 0.5 [dot(x_,r, x_,r) - dot(X_,r, X_,r)],
    //
    //     epsilon_ss = 0.5 [dot(x_,s, x_,s) - dot(X_,s, X_,s)],
    //
    //     gamma_rs = dot(x_,r, x_,s) - dot(X_,r, X_,s),
    //
    // and analogously for curvature and transverse shear.
    //
    // All reference quantities depend only on the undeformed geometry.
    // Consequently, they contribute no entries to B_nat.
    strain_nat(epsilon_rr) -= Precision(0.5) * point.X_rs.col(0).dot(point.X_rs.col(0));
    strain_nat(epsilon_ss) -= Precision(0.5) * point.X_rs.col(1).dot(point.X_rs.col(1));
    strain_nat(gamma_rs)   -= point.X_rs.col(0).dot(point.X_rs.col(1));
    strain_nat(kappa_rr)   -= point.X_rs.col(0).dot(point.D_rs.col(0));
    strain_nat(kappa_ss)   -= point.X_rs.col(1).dot(point.D_rs.col(1));
    strain_nat(kappa_rs)   -= point.X_rs.col(0).dot(point.D_rs.col(1)) + point.X_rs.col(1).dot(point.D_rs.col(0));
    strain_nat(gamma_r3)   -= point.X_rs.col(0).dot(point.D);
    strain_nat(gamma_s3)   -= point.X_rs.col(1).dot(point.D);
}

/**
 * Transforms generalized natural strain components into the pointwise local
 * orthonormal reference basis.
 *
 * The transformation depends exclusively on the undeformed curved-reference
 * geometry. It is therefore applied identically to values and B rows without
 * any additional nonlinear derivatives.
 *
 * @param point Pointwise curved-reference geometry.
 * @param strain Generalized strain vector transformed in place.
 * @param B Optional generalized B matrix transformed in place.
 */
template<Index N>
void FRTShell<N>::transform_strain_to_local(
    const ReferencePoint& point,
    Vec8&                 strain,
    Mat8x6N*              B
) const {
    const Precision t00 = point.invJ(0, 0);
    const Precision t01 = point.invJ(0, 1);
    const Precision t10 = point.invJ(1, 0);
    const Precision t11 = point.invJ(1, 1);

    StaticMatrix<3, 3> in_plane;
    in_plane << t00 * t00,               t01 * t01,               t00 * t01,
                t10 * t10,               t11 * t11,               t10 * t11,
                Precision(2) * t00 * t10, Precision(2) * t01 * t11, t00 * t11 + t01 * t10;

    // Transform the membrane and curvature engineering tensor blocks
    transform_in_plane_rows<N>(in_plane, 0, strain, B);
    transform_in_plane_rows<N>(in_plane, 3, strain, B);

    // Transform the covariant transverse-shear vector and the same B rows
    const Vec2 shear_nat = strain.template segment<2>(6);
    strain.template segment<2>(6) = point.invJ * shear_nat;

    if (B) {
        const Mat2x6N shear_B = B->template block<2, num_dofs>(6, 0);
        B->template block<2, num_dofs>(6, 0) = point.invJ * shear_B;
    }
}

/**
 * Prepares all pointwise quantities required for one shell evaluation.
 *
 * The caller controls which quantities are needed through the `with_*` flags.
 * Dependencies between these quantities are resolved at the beginning:
 *
 *     with_B          -> with_strain
 *     with_resultants -> with_strain, unless resultants are imported
 *
 * Temporary integration-point and tying-point data are stored in one
 * topology-specific `thread_local` workspace. The workspace belongs to the
 * executing thread, not to the element. Its vectors retain their allocated
 * capacity between calls, so repeated evaluations normally require no further
 * dynamic allocations.
 *
 * The returned `EvaluationData` object does not own this temporary storage. Its
 * spans point directly into the currently active parts of the thread-local
 * workspace and remain valid only until that workspace is reused by a later
 * evaluation on the same thread.
 *
 * Depending on the requested quantities, the function performs the following
 * steps:
 *
 *  1. Resize the required workspace buffers.
 *  2. Expose the active buffers through `EvaluationData`.
 *  3. Evaluate nodal SO(3) rotation derivatives when B or geometric terms are
 *     required.
 *  4. Initialize and scale the drilling stiffness.
 *  5. Evaluate the section tangent for linearized calls.
 *  6. Evaluate compatible strains and B matrices at tying points.
 *  7. Reconstruct MITC strains and B matrices at integration points.
 *  8. Import or evaluate generalized stress resultants.
 *
 * @param state Current or trial nodal shell state.
 * @param with_strain Request generalized strain values.
 * @param with_B Request first generalized strain derivatives.
 * @param with_G Request second rotation derivatives for the geometric tangent.
 * @param with_resultants Request generalized shell stress resultants.
 * @param ip_stress Optional previously evaluated integration-point resultants.
 * @param ip_start_idx First row of `ip_stress` belonging to this element.
 * @return Non-owning view of all quantities prepared in the active workspace.
 */
template<Index N>
typename FRTShell<N>::EvaluationData FRTShell<N>::init_evaluation(
    const CurrentState& state,
    bool                with_strain,
    bool                with_B,
    bool                with_G,
    bool                with_resultants,
    const Field*        ip_stress,
    int                 ip_start_idx
) const {
    // A B matrix is the first derivative of the generalized strain and
    // therefore cannot be evaluated without evaluating the strain kinematics.
    if (with_B) {
        with_strain = true;
    }

    // Constitutive resultants must be evaluated from the current strains unless
    // an already evaluated integration-point resultant field is supplied.
    if (with_resultants && ip_stress == nullptr) {
        with_strain = true;
    }

    const auto& ref = reference_data();

    const std::size_t num_ip    = ref.ip_points.size();
    const std::size_t num_tying = ref.tying_points.size();

    // One workspace exists for each executing thread and shell topology.
    //
    // The workspace is reused by all elements of the same topology processed
    // by that thread. Its vectors retain their capacity after resize(), so the
    // first sufficiently large evaluation establishes the temporary storage
    // required by later calls.
    thread_local EvaluationWorkspace workspace;

    // Resize only the buffers needed by the current evaluation.
    //
    // A size of zero marks a buffer as inactive while preserving its allocated
    // capacity for later calls.
    workspace.tying_strain_nat  .resize(with_strain     ? num_tying : 0);
    workspace.tying_B_nat       .resize(with_B          ? num_tying : 0);
    workspace.ip_strain         .resize(with_strain     ? num_ip : 0);
    workspace.ip_B              .resize(with_B          ? num_ip : 0);
    workspace.ip_resultants     .resize(with_resultants ? num_ip : 0);

    // The constitutive tangent is needed together with B. When resultants are
    // imported, no local material update is performed and no current tangent is
    // produced by compute_material_resultants().
    workspace.ip_tangent        .resize(with_B && ip_stress == nullptr ? num_ip : 0);
    workspace.ip_drill_stiffness.resize(with_B ? num_ip : 0);

    // Geometric tying weights depend on second kinematic derivatives and are
    // required only while assembling the geometric tangent.
    workspace.geometric_tying_weights.resize(with_G ? num_tying : 0);

    EvaluationData data;

    data.with_strain     = with_strain;
    data.with_B          = with_B;
    data.with_G          = with_G;
    data.with_resultants = with_resultants;
    data.state           = state;

    // Expose the active portions of the thread-local buffers.
    //
    // EvaluationData owns none of these arrays. The spans merely provide
    // convenient indexed access to the storage held by `workspace`.
    data.tying_strain_nat        = Span<Vec8>(workspace.tying_strain_nat);
    data.tying_B_nat             = Span<Mat8x6N>(workspace.tying_B_nat);

    data.ip_strain               = Span<Vec8>(workspace.ip_strain);
    data.ip_B                    = Span<Mat8x6N>(workspace.ip_B);
    data.ip_resultants           = Span<Vec8>(workspace.ip_resultants);
    data.ip_tangent              = Span<Mat8>(workspace.ip_tangent);

    data.ip_drill_stiffness      = Span<Precision>(workspace.ip_drill_stiffness);
    data.geometric_tying_weights = Span<Vec8>(workspace.geometric_tying_weights);

    // Evaluate nodal SO(3) derivatives once and reuse them at every tying and
    // integration point.
    //
    // with_B requires first derivatives
    //
    //     dR_i / d theta_i,a,
    //
    // while with_G additionally requires second derivatives
    //
    //     d²R_i / d theta_i,a d theta_i,b.
    //
    // The values are stored node-wise in the workspace because all pointwise
    // strain evaluations use the same nodal rotations.
    if (with_B || with_G) {
        for (Index node = 0; node < num_nodes; ++node) {
            const Vec3 theta = state.theta.row(node).transpose();

            if (with_G) {
                math::so3::rotation_matrix_second_derivatives(
                    theta,
                    workspace.rotations[node].value,
                    workspace.rotations[node].d1,
                    workspace.rotations[node].d2
                );
            } else {
                math::so3::rotation_matrix_first_derivatives(
                    theta,
                    workspace.rotations[node].value,
                    workspace.rotations[node].d1
                );
            }
        }

        // Store a non-owning pointer to the complete nodal rotation cache.
        data.rotations = &workspace.rotations;
    }

    // Initialize the configuration-independent drilling modulus once for every
    // element.
    //
    // The underlying modulus is derived from the zero-strain in-plane shear
    // tangent
    //
    //     A66 = H0(gamma_12, gamma_12).
    //
    // Only the positive magnitude is retained. Keeping this base value
    // independent of the current deformation ensures that drilling force and
    // drilling tangent derive from the same fixed quadratic potential.
    if (with_B) {
        std::call_once(ref.drill_stiffness_once, [&]() {
            using Component = ShellGeneralizedStrain::Component;

            constexpr Index gamma_12 =
                static_cast<Index>(Component::GammaXY);

            for (ReferencePoint& point : reference_data_->ip_points) {
                ShellGeneralizedStrain zero_strain(Vec8::Zero());
                ShellStressResultants  zero_resultants;
                Mat8                   H0;

                Mat3 basis = point.basis;

                shell_section()->evaluate(
                    reference_position(point.r, point.s),
                    basis,
                    zero_strain,
                    false,
                    zero_resultants,
                    H0
                );

                point.drill_stiffness_base =
                    drill_scale * std::abs(H0(gamma_12, gamma_12));
            }
        });

        // Apply the current element-topology scale to the cached base modulus.
        //
        // The base modulus is stored in the reference data, whereas the scaled
        // pointwise values belong to the current evaluation workspace.
        const Precision stiffness_scale =
            std::abs(topology_stiffness_scale());

        for (Index ip = 0; ip < static_cast<Index>(num_ip); ++ip) {
            const std::size_t id = static_cast<std::size_t>(ip);

            data.ip_drill_stiffness[id] =
                stiffness_scale * ref.ip_points[id].drill_stiffness_base;
        }
    }

    // Linearized evaluations require the zero-strain shell-section tangent at
    // every integration point.
    //
    // Nonlinear evaluations with resultants obtain their current constitutive
    // tangent later from compute_material_resultants(). Therefore this explicit
    // tangent evaluation is needed only when B is requested without resultants.
    if (with_B && !with_resultants) {
        for (Index ip = 0; ip < static_cast<Index>(num_ip); ++ip) {
            const std::size_t    id    = static_cast<std::size_t>(ip);
            const ReferencePoint& point = ref.ip_points[id];

            data.ip_tangent[id] = resultant_stiffness(point.r, point.s);
        }
    }

    if (with_strain) {
        // First evaluate the compatible natural strains at all tying points.
        //
        // MITC interpolation at the integration points depends on these tying
        // values, so this loop must be completed before the integration-point
        // loop begins.
        for (Index tying = 0;
             tying < static_cast<Index>(num_tying);
             ++tying) {
            const std::size_t id = static_cast<std::size_t>(tying);

            compute_natural_strain(
                data,
                ref.tying_points[id],
                data.tying_strain_nat[id],
                with_B ? &data.tying_B_nat[id] : nullptr
            );
        }

        // Evaluate the complete strain field at every integration point:
        //
        //  1. Compute the compatible strain and its optional B matrix in
        //     natural coordinates.
        //  2. Replace the MITC-controlled components by values reconstructed
        //     from the tying points.
        //  3. Transform all generalized strain components and B rows from the
        //     natural coordinates into the local material basis.
        for (Index ip = 0; ip < static_cast<Index>(num_ip); ++ip) {
            const std::size_t     id    = static_cast<std::size_t>(ip);
            const ReferencePoint& point = ref.ip_points[id];

            Vec8&    strain = data.ip_strain[id];
            Mat8x6N* B      = with_B ? &data.ip_B[id] : nullptr;

            compute_natural_strain(data, point, strain, B);
            apply_mitc_natural(data, point, strain, B);
            transform_strain_to_local(point, strain, B);
        }
    }

    if (with_resultants) {
        if (ip_stress) {
            // Reuse generalized resultants supplied by the caller. No local
            // constitutive section update is performed in this path.
            load_ip_resultants(data, *ip_stress, ip_start_idx);
        } else {
            // Evaluate generalized stress resultants and the current
            // constitutive tangent from the integration-point strains.
            compute_material_resultants(data);
        }
    }

    return data;
}

} // namespace fem::model
