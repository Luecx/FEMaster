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
 * Adds one scaled scalar product of two current midsurface derivatives to a
 * membrane strain component and optionally to its B row.
 *
 * The current vector fields are
 *
 *     x_a = sum_i a_i x_i,
 *     x_b = sum_i b_i x_i.
 *
 * Their scalar product contributes `scale * x_a . x_b`. Because nodal
 * positions depend linearly on the translational degrees of freedom, the
 * derivative is inserted directly as
 *
 *     d(x_a . x_b)/du_i = a_i x_b + b_i x_a.
 *
 * No explicit `dx/du` identity matrices are constructed or stored.
 *
 * @tparam N Number of shell nodes.
 * @param positions Current nodal midsurface positions.
 * @param a_coefficients Interpolation coefficients of the first vector field.
 * @param b_coefficients Interpolation coefficients of the second vector field.
 * @param scale Scalar multiplier of the complete contribution.
 * @param strain_id Generalized natural strain row receiving the contribution.
 * @param strain Generalized natural strain vector to update.
 * @param B Optional generalized strain-displacement matrix to update.
 */
template<Index N>
void add_xx_strain_and_B(
    const typename FRTShell<N>::MatN3& positions,
    const typename FRTShell<N>::VecN&  a_coefficients,
    const typename FRTShell<N>::VecN&  b_coefficients,
    Precision                          scale,
    Index                              strain_id,
    typename FRTShell<N>::Vec8&        strain,
    typename FRTShell<N>::Mat8x6N*     B
) {
    using Shell = FRTShell<N>;

    Vec3 a_value = Vec3::Zero();
    Vec3 b_value = Vec3::Zero();

    // Interpolate the two current vector fields from the nodal positions
    for (Index node = 0; node < N; ++node) {
        const Vec3 x_i = positions.row(node).transpose();
        a_value += a_coefficients(node) * x_i;
        b_value += b_coefficients(node) * x_i;
    }

    // Add the scalar kinematic value to the requested strain component
    strain(strain_id) += scale * a_value.dot(b_value);

    if (!B) {
        return;
    }

    // Insert the analytical translational derivative directly into the B row
    for (Index node = 0; node < N; ++node) {
        const Index base = Shell::dofs_per_node * node;
        const Vec3 derivative = scale
                              * (a_coefficients(node) * b_value
                               + b_coefficients(node) * a_value);

        B->template block<1, 3>(strain_id, base) += derivative.transpose();
    }
}

/**
 * Adds one scaled scalar product between a current midsurface derivative and a
 * current director field to a generalized strain component and optional B row.
 *
 * The two vector fields are
 *
 *     x_a = sum_i a_i x_i,
 *     d_b = sum_i b_i R_i d0_i.
 *
 * Translational derivatives are inserted analytically. Rotational derivatives
 * use the compact nodal matrices `dR_i/dtheta_ia` retained by the active
 * thread-local workspace:
 *
 *     d(d_b)/dtheta_ia = b_i (dR_i/dtheta_ia) d0_i.
 *
 * @tparam N Number of shell nodes.
 * @param data Active element evaluation data.
 * @param reference_directors Nodal undeformed directors acted on by R.
 * @param a_coefficients Interpolation coefficients of the position field.
 * @param b_coefficients Interpolation coefficients of the director field.
 * @param scale Scalar multiplier of the complete contribution.
 * @param strain_id Generalized natural strain row receiving the contribution.
 * @param strain Generalized natural strain vector to update.
 * @param B Optional generalized strain-displacement matrix to update.
 */
template<Index N>
void add_xd_strain_and_B(
    const typename FRTShell<N>::EvaluationData& data,
    const typename FRTShell<N>::MatN3&          reference_directors,
    const typename FRTShell<N>::VecN&           a_coefficients,
    const typename FRTShell<N>::VecN&           b_coefficients,
    Precision                                   scale,
    Index                                       strain_id,
    typename FRTShell<N>::Vec8&                 strain,
    typename FRTShell<N>::Mat8x6N*              B
) {
    using Shell = FRTShell<N>;

    Vec3 x_value = Vec3::Zero();
    Vec3 d_value = Vec3::Zero();

    // Interpolate the current position derivative and director field once
    for (Index node = 0; node < N; ++node) {
        x_value += a_coefficients(node) * data.state.x.row(node).transpose();
        d_value += b_coefficients(node) * data.state.d.row(node).transpose();
    }

    // Add the current scalar product to the requested strain component
    strain(strain_id) += scale * x_value.dot(d_value);

    if (!B) {
        return;
    }

    logging::error(data.rotations != nullptr,
                   "FRTShell: B evaluation requires nodal rotation derivatives");

    const auto& rotations = *data.rotations;

    // Differentiate the position field with respect to nodal translations
    for (Index node = 0; node < N; ++node) {
        const Index base = Shell::dofs_per_node * node;
        const Vec3 derivative = scale * a_coefficients(node) * d_value;
        B->template block<1, 3>(strain_id, base) += derivative.transpose();
    }

    // Differentiate the director field with respect to the three local nodal
    // axis-angle coordinates using the persistent nodal reference director.
    for (Index node = 0; node < N; ++node) {
        const Index     rot_base    = Shell::dofs_per_node * node + 3;
        const Precision coefficient = b_coefficients(node);

        if (coefficient == Precision(0)) {
            continue;
        }

        const Vec3 d0 = reference_directors.row(node).transpose();

        for (Index component = 0; component < 3; ++component) {
            const Vec3 director_derivative = rotations[node].d1[component] * d0;
            (*B)(strain_id, rot_base + component) +=
                scale * coefficient * director_derivative.dot(x_value);
        }
    }
}

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
 * Membrane strains are Green-Lagrange midsurface strains. Bending strains
 * compare current position/director derivatives with their reference values,
 * and transverse-shear strains compare current tangents with the interpolated
 * current director. Engineering shear conventions are used for the `rs`,
 * curvature and transverse-shear entries.
 *
 * When `B_nat` is supplied, all first derivatives are assembled analytically.
 * Trivial translational identities are inserted directly, while director
 * derivatives use only compact nodal SO(3) matrices.
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

    // Copy the pointwise interpolation coefficients once because each group of
    // shell strains reuses the same natural derivatives and shape values
    const VecN dshape_r = point.dshape_rs.col(0);
    const VecN dshape_s = point.dshape_rs.col(1);
    const VecN shape    = point.shape;
    const MatN3& d0     = reference_data().d0;

    // Construct the three membrane metric changes
    add_xx_strain_and_B<N>(
        data.state.x,
        dshape_r,
        dshape_r,
        Precision(0.5),
        epsilon_rr,
        strain_nat,
        B_nat
    );
    add_xx_strain_and_B<N>(
        data.state.x,
        dshape_s,
        dshape_s,
        Precision(0.5),
        epsilon_ss,
        strain_nat,
        B_nat
    );
    add_xx_strain_and_B<N>(
        data.state.x,
        dshape_r,
        dshape_s,
        Precision(1),
        gamma_rs,
        strain_nat,
        B_nat
    );

    // Construct the three changes of curvature from current tangent/director
    // derivative products
    add_xd_strain_and_B<N>(data, d0, dshape_r, dshape_r, Precision(1), kappa_rr, strain_nat, B_nat);
    add_xd_strain_and_B<N>(data, d0, dshape_s, dshape_s, Precision(1), kappa_ss, strain_nat, B_nat);
    add_xd_strain_and_B<N>(data, d0, dshape_r, dshape_s, Precision(1), kappa_rs, strain_nat, B_nat);
    add_xd_strain_and_B<N>(data, d0, dshape_s, dshape_r, Precision(1), kappa_rs, strain_nat, B_nat);

    // Construct the two transverse-shear measures from current tangents and the
    // interpolated current director
    add_xd_strain_and_B<N>(data, d0, dshape_r, shape, Precision(1), gamma_r3, strain_nat, B_nat);
    add_xd_strain_and_B<N>(data, d0, dshape_s, shape, Precision(1), gamma_s3, strain_nat, B_nat);

    // Subtract the complete reference metric, curvature and shear values so the
    // undeformed curved shell has exactly zero generalized strain
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
 * Prepares all quantities requested by one shell evaluation.
 *
 * The function owns no dynamic temporary data. A `thread_local` workspace is
 * resized to the required topology-specific point counts and retains its
 * capacity for subsequent elements on the same worker thread. The returned
 * `EvaluationData` contains non-owning spans into this workspace.
 *
 * First and second SO(3) derivatives are evaluated once per node and reused at
 * all tying and integration points. The cached zero-strain in-plane shear
 * modulus A66 defines a configuration-independent drilling coefficient at every
 * integration point, ensuring that force and tangent derive from one fixed
 * quadratic potential even for nonlinear material sections.
 *
 * @param state Current or trial nodal shell state.
 * @param with_strain Request generalized strain values.
 * @param with_B Request first generalized strain derivatives.
 * @param with_G Request second nodal SO(3) derivatives.
 * @param with_resultants Request generalized shell resultants.
 * @param ip_stress Optional previously evaluated resultant field.
 * @param ip_start_idx First row of the supplied resultant field.
 * @return Non-owning view into the active thread-local workspace.
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
    if (with_B) {
        with_strain = true;
    }

    if (with_resultants && ip_stress == nullptr) {
        with_strain = true;
    }

    const auto& ref = reference_data();
    const std::size_t num_ip    = ref.ip_points.size();
    const std::size_t num_tying = ref.tying_points.size();

    // One reusable workspace exists for every executing thread and shell
    // topology. No element owns or deallocates this temporary storage.
    thread_local EvaluationWorkspace workspace;

    // Resize only the arrays required by this caller. Standard vector resize
    // preserves capacity, so repeated evaluations do not allocate after the
    // thread has encountered its largest element request.
    workspace.tying_strain_nat.resize(with_strain ? num_tying : 0);
    workspace.tying_B_nat.resize(with_B ? num_tying : 0);
    workspace.ip_strain.resize(with_strain ? num_ip : 0);
    workspace.ip_B.resize(with_B ? num_ip : 0);
    workspace.ip_resultants.resize(with_resultants ? num_ip : 0);
    workspace.ip_tangent.resize(with_B && ip_stress == nullptr ? num_ip : 0);
    workspace.ip_drill_stiffness.resize(with_B ? num_ip : 0);
    workspace.geometric_tying_weights.resize(with_G ? num_tying : 0);

    EvaluationData data;
    data.with_strain     = with_strain;
    data.with_B          = with_B;
    data.with_G          = with_G;
    data.with_resultants = with_resultants;
    data.state           = state;

    // Expose exactly the active portions of the thread-local buffers
    data.tying_strain_nat = Span<Vec8>(workspace.tying_strain_nat);
    data.tying_B_nat = Span<Mat8x6N>(workspace.tying_B_nat);
    data.ip_strain = Span<Vec8>(workspace.ip_strain);
    data.ip_B = Span<Mat8x6N>(workspace.ip_B);
    data.ip_resultants = Span<Vec8>(workspace.ip_resultants);
    data.ip_tangent = Span<Mat8>(workspace.ip_tangent);
    data.ip_drill_stiffness = Span<Precision>(workspace.ip_drill_stiffness);
    data.geometric_tying_weights = Span<Vec8>(workspace.geometric_tying_weights);

    // Evaluate each nodal SO(3) rotation and its requested derivatives once.
    // The complete local 3 x 3 matrices are retained because both the physical
    // director and objective drilling tangent vectors act on the same rotation.
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

        data.rotations = &workspace.rotations;
    }

    // Initialize the unscaled positive drilling modulus lazily and exactly once
    // for this element. This keeps temporary-storage lifecycle independent of
    // step_begin()/step_end() while avoiding repeated zero-strain section calls.
    if (with_B) {
        std::call_once(ref.drill_stiffness_once, [&]() {
            using Component = ShellGeneralizedStrain::Component;
            constexpr Index gamma_12 = static_cast<Index>(Component::GammaXY);

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

        // Apply only the current topology scale during each evaluation. The
        // underlying modulus remains independent of the current deformation,
        // so force and tangent derive from one fixed quadratic potential.
        const Precision stiffness_scale = std::abs(topology_stiffness_scale());

        for (Index ip = 0; ip < static_cast<Index>(num_ip); ++ip) {
            data.ip_drill_stiffness[static_cast<std::size_t>(ip)] =
                stiffness_scale
                * ref.ip_points[static_cast<std::size_t>(ip)].drill_stiffness_base;
        }
    }

    // Linearized callers without a constitutive material update require the
    // zero-strain section tangent at every integration point. Nonlinear callers
    // receive their current tangent from compute_material_resultants().
    if (with_B && !with_resultants) {
        for (Index ip = 0; ip < static_cast<Index>(num_ip); ++ip) {
            const ReferencePoint& point = ref.ip_points[static_cast<std::size_t>(ip)];
            data.ip_tangent[static_cast<std::size_t>(ip)] =
                resultant_stiffness(point.r, point.s);
        }
    }

    // Tying values must be prepared before the integration-point MITC fields
    // can be reconstructed. The two loops remain here because they form one
    // evaluation sequence and are not reused independently anywhere else.
    if (with_strain) {
        for (Index tying = 0; tying < static_cast<Index>(ref.tying_points.size()); ++tying) {
            const std::size_t id = static_cast<std::size_t>(tying);
            compute_natural_strain(
                data,
                ref.tying_points[id],
                data.tying_strain_nat[id],
                data.with_B ? &data.tying_B_nat[id] : nullptr
            );
        }

        for (Index ip = 0; ip < static_cast<Index>(ref.ip_points.size()); ++ip) {
            const std::size_t id = static_cast<std::size_t>(ip);
            const ReferencePoint& point = ref.ip_points[id];
            Vec8& strain = data.ip_strain[id];
            Mat8x6N* B = data.with_B ? &data.ip_B[id] : nullptr;

            // Evaluate the compatible natural field, replace the selected MITC
            // components and finally transform to the local material basis
            compute_natural_strain(data, point, strain, B);
            apply_mitc_natural(data, point, strain, B);
            transform_strain_to_local(point, strain, B);
        }
    }

    // Either import an existing resultant field or perform the nonlinear
    // constitutive shell-section update at all integration points
    if (with_resultants) {
        if (ip_stress) {
            load_ip_resultants(data, *ip_stress, ip_start_idx);
        } else {
            compute_material_resultants(data);
        }
    }

    return data;
}

} // namespace fem::model
