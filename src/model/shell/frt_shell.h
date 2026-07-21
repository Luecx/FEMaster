/**
 * @file frt_shell.h
 * @brief Declares the common finite-rotation Total-Lagrangian shell element.
 *
 * The class implements the topology-independent part of geometrically
 * nonlinear Reissner-Mindlin shell elements with three translational and three
 * rotational degrees of freedom per node. Concrete S3, S4, S6 and S8 elements
 * provide their shape functions, quadrature rules, surface representation and
 * MITC assumed-strain interpolation.
 *
 * The formulation integrates over the actual isoparametric reference
 * midsurface. A local orthonormal reference basis is constructed independently
 * at every integration, tying and output point. No planar substitute geometry
 * is introduced for curved or warped shell elements.
 *
 * The generalized shell strain vector is ordered as
 *
 *     [eps11, eps22, gamma12, kappa11, kappa22, kappa12,
 *      gamma13, gamma23].
 *
 * The consistent nonlinear tangent is assembled without storing one dense
 * element Hessian per generalized strain component. The material part uses
 *
 *     K_mat = integral_A0 B^T H B dA0,
 *
 * while the geometric part contracts the generalized resultants directly with
 * the structured translational and rotational strain Hessian blocks.
 *
 * Drilling rotation is stabilized through an objective in-plane spin measure.
 * The stabilization vanishes exactly for arbitrary finite rigid-body motion
 * and is derived from one quadratic potential, so its force and tangent remain
 * mutually consistent.
 *
 * Temporary evaluation arrays are retained in one `thread_local` workspace per
 * shell topology and worker thread. Element objects therefore carry only their
 * persistent reference geometry, while repeated element evaluations avoid
 * dynamic allocation after the workspace has reached its maximum required
 * capacity.
 *
 * @see FRTShellS3
 * @see FRTShellS4
 * @see FRTShellS6
 * @see FRTShellS8
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#pragma once

#include "shell.h"

#include "../../math/quadrature.h"

#include <array>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace fem::model {

/**
 * @brief Common finite-rotation Total-Lagrangian Reissner-Mindlin shell.
 *
 * `FRTShell<N>` contains the nonlinear nodal rotation kinematics, exact first
 * and second SO(3) derivatives, pointwise curved-reference geometry,
 * constitutive shell-section evaluation, consistent element assembly,
 * objective drilling stabilization, mass integration and result recovery
 * shared by all supported shell topologies.
 *
 * Every node carries the global degrees of freedom
 * `[ux, uy, uz, rx, ry, rz]`. The rotational coordinates form one total
 * axis-angle vector. The current nodal rotation is evaluated through the SO(3)
 * exponential map. Only the compact local `3 x 3` rotation derivatives are
 * retained; trivial translational identities such as `dx_i/du_i = I` are
 * inserted analytically where they are required.
 *
 * Concrete elements must provide their interpolation, quadrature and MITC
 * operators. The common class applies the same assumed-strain interpolation to
 * values and B rows and uses the transposed operator for direct geometric-
 * tangent assembly.
 *
 * @tparam N Number of shell midsurface nodes.
 */
template<Index N>
struct FRTShell : ShellElement<N> {
    using Base = ShellElement<N>;

    static_assert(N == 3 || N == 4 || N == 6 || N == 8,
                  "FRTShell supports only S3, S4, S6 and S8 topologies");

    // Number of midsurface interpolation nodes belonging to this topology.
    static constexpr Index num_nodes = N;

    // Number of global translational and rotational degrees of freedom stored
    // at every shell node.
    static constexpr Index dofs_per_node = 6;

    // Total number of element degrees of freedom in nodal ordering.
    static constexpr Index num_dofs = num_nodes * dofs_per_node;

    // Number of generalized membrane, bending and transverse-shear components.
    static constexpr Index num_strains = 8;

    // Relative stiffness of the objective drilling penalty with respect to the
    // zero-strain in-plane shear resultant stiffness A66. The coefficient is
    // deliberately small because the term only removes the rank deficiency of
    // the director-only Reissner-Mindlin kinematics.
    static constexpr Precision drill_scale = Precision(1e-5);

    // Relative artificial drilling inertia with respect to the physical rotary
    // inertia rho*h^3/12. The scaling is dimensionless and therefore preserves
    // the correct rotational mass dimension kg*m^2.
    static constexpr Precision drill_inertia_scale = Precision(1e-6);

    // Fixed-size generalized strain or resultant vector.
    using Vec8 = StaticVector<num_strains>;

    // Fixed-size element vector in nodal degree-of-freedom ordering.
    using Vec6N = StaticVector<num_dofs>;

    // Fixed-size generalized shell-section tangent matrix.
    using Mat8 = StaticMatrix<num_strains, num_strains>;

    // Fixed-size element stiffness or mass matrix.
    using Mat6N = StaticMatrix<num_dofs, num_dofs>;

    // Fixed-size generalized strain-displacement matrix.
    using Mat8x6N = StaticMatrix<num_strains, num_dofs>;

    // Two-row matrix used for transverse-shear transformations.
    using Mat2x6N = StaticMatrix<2, num_dofs>;

    // Three-row matrix used for vector-valued element derivatives.
    using Mat3x6N = StaticMatrix<3, num_dofs>;

    // Nodal three-dimensional vector field stored row-wise.
    using MatN3 = StaticMatrix<num_nodes, 3>;

    // Nodal natural-coordinate field stored row-wise.
    using MatN2 = StaticMatrix<num_nodes, 2>;

    // Nodal six-component position and rotation field stored row-wise.
    using MatN6 = StaticMatrix<num_nodes, 6>;

    // Fixed-size vector containing one scalar coefficient per element node.
    using VecN = StaticVector<num_nodes>;

    /**
     * @brief Compact value and derivatives of one nodal SO(3) rotation.
     *
     * The structure stores the complete nodal rotation matrix because the shell
     * director and the two drilling tangent vectors require the same rotation
     * but act on different reference vectors. The derivative convention is
     *
     *     d1[a]    = dR / d(theta_a),
     *     d2[a][b] = d²R / (d(theta_a) d(theta_b)).
     *
     * Every matrix is only `3 x 3`. No element-sized sparse-looking derivative
     * matrix is materialized. Second derivatives are written only when the
     * caller requests the nonlinear geometric tangent and are never read
     * otherwise.
     */
    struct RotationDerivatives {
        // Current nodal rotation matrix R(theta).
        Mat3 value;

        // First derivatives of R with respect to the three nodal axis-angle
        // coordinates.
        std::array<Mat3, 3> d1;

        // Second derivatives of R with respect to pairs of nodal axis-angle
        // coordinates.
        std::array<std::array<Mat3, 3>, 3> d2;
    };

    /**
     * @brief Current nodal positions, directors and total rotation vectors.
     *
     * `x` contains the current midsurface positions, `d` contains the rotated
     * nodal directors and `theta` contains the global total axis-angle vectors.
     * The state is intentionally compact and does not store derivative
     * identities or element-sized auxiliary matrices.
     */
    struct CurrentState {
        // Current global midsurface position at every element node.
        MatN3 x = MatN3::Zero();

        // Current global director at every element node.
        MatN3 d = MatN3::Zero();

        // Current total global axis-angle vector at every element node.
        MatN3 theta = MatN3::Zero();
    };

    /**
     * @brief Reference geometry evaluated at one natural shell point.
     *
     * The point stores the shape functions, actual curved reference tangents,
     * surface Jacobian, local orthonormal basis and interpolated reference
     * director field. `A` maps local tangent derivatives `(a,b)` to natural
     * derivatives `(r,s)`:
     *
     *     [N_,r, N_,s]^T = A [N_,a, N_,b]^T.
     */
    struct ReferencePoint {
        // First natural coordinate of the point.
        Precision r = Precision(0);

        // Second natural coordinate of the point.
        Precision s = Precision(0);

        // Numerical quadrature weight; zero for tying and arbitrary output
        // points.
        Precision w = Precision(0);

        // Physical reference-surface Jacobian |X_,r x X_,s|.
        Precision detJ = Precision(0);

        // Shape-function values at the point.
        VecN shape = VecN::Zero();

        // Natural shape-function derivatives [N_,r, N_,s].
        MatN2 dshape_rs = MatN2::Zero();

        // Shape-function derivatives along the first orthonormal reference
        // tangent.
        VecN dshape_a = VecN::Zero();

        // Shape-function derivatives along the second orthonormal reference
        // tangent.
        VecN dshape_b = VecN::Zero();

        // First natural reference midsurface tangent X_,r.
        Vec3 X_r = Vec3::Zero();

        // Second natural reference midsurface tangent X_,s.
        Vec3 X_s = Vec3::Zero();

        // First orthonormal reference tangent X_,a = e1.
        Vec3 X_a = Vec3::Zero();

        // Second orthonormal reference tangent X_,b = e2.
        Vec3 X_b = Vec3::Zero();

        // Interpolated nodal reference director field.
        Vec3 D = Vec3::Zero();

        // Natural derivative D_,r of the reference director field.
        Vec3 D_r = Vec3::Zero();

        // Natural derivative D_,s of the reference director field.
        Vec3 D_s = Vec3::Zero();

        // Local tangent derivative D_,a of the reference director field.
        Vec3 D_a = Vec3::Zero();

        // Local tangent derivative D_,b of the reference director field.
        Vec3 D_b = Vec3::Zero();

        // First pointwise orthonormal reference tangent.
        Vec3 e1 = Vec3::Zero();

        // Second pointwise orthonormal reference tangent.
        Vec3 e2 = Vec3::Zero();

        // Pointwise reference surface normal.
        Vec3 e3 = Vec3::Zero();

        // Mapping from local tangent derivatives to natural derivatives.
        Mat2 A = Mat2::Zero();

        // Mapping from natural derivatives to local tangent derivatives.
        Mat2 invA = Mat2::Zero();

        // Unscaled objective drilling stiffness per unit reference area. The
        // value equals drill_scale times the zero-strain in-plane shear
        // resultant stiffness A66 and is cached once with the reference point.
        Precision drill_stiffness_base = Precision(0);
    };

    /**
     * @brief Persistent reference configuration of one shell element.
     *
     * The data is initialized once for the active analysis state and remains
     * unchanged during all Newton iterations. Integration and tying points use
     * the actual isoparametric reference surface, including curved S6/S8 and
     * warped S4 geometry.
     */
    struct ReferenceData {
        // Global reference midsurface coordinates of the element nodes.
        MatN3 X = MatN3::Zero();

        // Global reference directors of the element nodes.
        MatN3 d0 = MatN3::Zero();

        // Cached reference geometry at all numerical integration points.
        std::vector<ReferencePoint> ip_points;

        // Cached reference geometry at all topology-specific MITC tying points.
        std::vector<ReferencePoint> tying_points;

        // Reference midsurface area obtained from the selected quadrature rule.
        Precision area = Precision(0);

        // One-time guard for the configuration-independent drilling moduli.
        // The moduli are initialized lazily on the first evaluation that
        // requests a B matrix, independently of step_begin()/step_end().
        mutable std::once_flag drill_stiffness_once;
    };

    /**
     * @brief Reusable temporary storage owned independently by every thread.
     *
     * One workspace exists per template topology and executing thread. The
     * vectors retain their capacity between element calls, eliminating repeated
     * heap allocation while avoiding any memory growth proportional to the
     * number of shell elements. The workspace is intentionally non-reentrant:
     * one thread may have only one active `FRTShell<N>` evaluation at a time.
     */
    struct EvaluationWorkspace {
        // Compact nodal rotations and their local SO(3) derivatives.
        std::array<RotationDerivatives, num_nodes> rotations;

        // Compatible natural generalized strains at all MITC tying points.
        std::vector<Vec8> tying_strain_nat;

        // Compatible natural strain-displacement matrices at all tying points.
        std::vector<Mat8x6N> tying_B_nat;

        // Local generalized strains at all integration points.
        std::vector<Vec8> ip_strain;

        // Local generalized strain-displacement matrices at all integration
        // points.
        std::vector<Mat8x6N> ip_B;

        // Local generalized shell resultants at all integration points.
        std::vector<Vec8> ip_resultants;

        // Local generalized shell-section tangents at all integration points.
        std::vector<Mat8> ip_tangent;

        // Configuration-independent drilling stiffness per unit reference area
        // at every integration point.
        std::vector<Precision> ip_drill_stiffness;

        // Reusable MITC pull-back weights for all tying points. The buffer is
        // reset and reused for every integration point during geometric-tangent
        // assembly.
        std::vector<Vec8> geometric_tying_weights;
    };

    /**
     * @brief Minimal C++17 non-owning contiguous view.
     *
     * This mirrors the small subset of `std::span` needed by the shell
     * workspaces while keeping the project on its C++17 baseline.
     */
    template<typename T>
    struct Span {
        T*          data_ = nullptr;
        std::size_t size_ = 0;

        Span() = default;

        explicit Span(std::vector<T>& values)
            : data_(values.data()),
              size_(values.size()) {}

        T& operator[](std::size_t index) const {
            return data_[index];
        }

        T* begin() const {
            return data_;
        }

        T* end() const {
            return data_ + size_;
        }

        std::size_t size() const {
            return size_;
        }

        bool empty() const {
            return size_ == 0;
        }
    };

    /**
     * @brief Non-owning view of all quantities required by one element call.
     *
     * The spans refer to the current thread's `EvaluationWorkspace`. They are
     * valid only until the next shell evaluation of the same topology starts on
     * that thread. No pointer or span may be retained by an element, section or
     * solver after the active element call returns.
     */
    struct EvaluationData {
        // Whether generalized strain values were requested.
        bool with_strain = false;

        // Whether first strain derivatives were requested.
        bool with_B = false;

        // Whether second nodal rotation derivatives were requested for the
        // geometric tangent.
        bool with_G = false;

        // Whether generalized shell resultants were requested.
        bool with_resultants = false;

        // Compact current nodal configuration used by all kinematic routines.
        CurrentState state;

        // Nodal SO(3) values and derivatives retained in the thread-local
        // workspace. The pointer is null when no rotational derivatives are
        // required.
        const std::array<RotationDerivatives, num_nodes>* rotations = nullptr;

        // Compatible natural generalized strains at MITC tying points.
        Span<Vec8> tying_strain_nat;

        // Compatible natural B matrices at MITC tying points.
        Span<Mat8x6N> tying_B_nat;

        // Local generalized strains at numerical integration points.
        Span<Vec8> ip_strain;

        // Local generalized B matrices at numerical integration points.
        Span<Mat8x6N> ip_B;

        // Local generalized resultants at numerical integration points.
        Span<Vec8> ip_resultants;

        // Local generalized section tangents at numerical integration points.
        Span<Mat8> ip_tangent;

        // Objective drilling penalty stiffness per unit reference area at every
        // integration point.
        Span<Precision> ip_drill_stiffness;

        // Reusable tying-point buffer used by the transposed MITC interpolation
        // during geometric-tangent assembly.
        Span<Vec8> geometric_tying_weights;
    };

    // Persistent curved-reference geometry belonging only to this element.
    std::unique_ptr<ReferenceData> reference_data_;

    // Lifecycle and element construction.
    // These functions own the persistent per-element state. The constructor
    // stores the element id/connectivity through the shell base class, while
    // step_begin() builds the reference geometry used by all later integration
    // calls. step_end() releases only that persistent reference cache; reusable
    // thread-local evaluation workspaces are intentionally independent of the
    // element lifecycle.
    FRTShell(ID id, const std::array<ID, N>& nodes);
    ~FRTShell() override = default;
    void step_begin() override;
    void step_end() override;

    // Topology-specific hooks implemented by the concrete MITC elements.
    // The common FRT shell code knows the nonlinear rotation kinematics,
    // reference mapping, material calls and assembly loops. The concrete S3,
    // S4, S6 and S8 classes provide only the data that depends on nodal
    // topology: shape functions, natural node coordinates, quadrature point
    // locations exposed for output, surface extraction, and the forward and
    // transposed MITC assumed-strain operators.
    virtual VecN shape_function(Precision r, Precision s) const = 0;
    virtual MatN2 shape_derivative(Precision r, Precision s) const = 0;
    virtual MatN2 node_coords_natural() const = 0;
    virtual std::vector<Vec2> tying_point_coordinates() const = 0;
    virtual void apply_mitc_natural(
        const EvaluationData& data,
        const ReferencePoint& point,
        Vec8&                 strain_nat,
        Mat8x6N*              B_nat
    ) const = 0;
    virtual void pull_back_mitc_resultants(
        const ReferencePoint& point,
        const Vec8&           assumed_weights,
        Vec8&                 compatible_weights,
        Span<Vec8>            tying_weights
    ) const = 0;

    // Reference geometry helpers.
    // Everything in this block is based on the undeformed isoparametric shell
    // midsurface. The helpers initialize nodal reference coordinates and
    // directors, construct pointwise curved reference data, and retrieve cached
    // integration or tying points. These functions deliberately stay separate
    // from current-state kinematics so reference data can be reused across all
    // Newton iterations of the active step.
    const ReferenceData& reference_data() const;
    MatN3 init_reference_node_coords() const;
    MatN3 init_reference_directors(const MatN3& X) const;
    ReferencePoint make_reference_point(Precision r, Precision s, Precision w) const;
    Vec3 reference_position(Precision r, Precision s) const;
    Mat3 reference_basis_global(Precision r, Precision s) const;
    const ReferencePoint* cached_reference_point(Precision r, Precision s) const;

    // Current configuration and displacement extraction.
    // These functions gather element-local state from model fields. The current
    // POSITION field contains translational coordinates and total axis-angle
    // rotations, while displacement-based helpers build trial states from the
    // reference geometry plus a supplied displacement field. All returned data
    // is compact nodal data; derivative workspaces are filled later only when an
    // evaluation explicitly requests them.
    MatN6 node_coords_current_6() const;
    CurrentState current_state() const;
    CurrentState reference_state() const;
    CurrentState current_state_from_displacement(const Field& displacement) const;
    Vec6N element_displacement_vector(const Field& displacement) const;

    // Kinematic strain evaluation.
    // This block builds the compatible natural strain field and maps it into
    // the local orthonormal reference basis used by shell sections. The same
    // path is used for values and B rows, so MITC replacement and basis
    // transformation remain consistent for linear stiffness, nonlinear force,
    // geometric tangent and result recovery.
    void compute_natural_strain(
        const EvaluationData& data,
        const ReferencePoint& point,
        Vec8&                 strain_nat,
        Mat8x6N*              B_nat = nullptr
    ) const;
    void transform_strain_to_local(
        const ReferencePoint& point,
        Vec8&                 strain,
        Mat8x6N*              B = nullptr
    ) const;

    // Evaluation workspace preparation and shell-section response.
    // init_evaluation() is the central setup step for one element call. It
    // resizes the thread-local buffers, exposes active spans through
    // EvaluationData, evaluates requested SO(3) derivatives, prepares tying and
    // integration-point strains, and either imports stored resultants or asks
    // the shell section for fresh nonlinear resultants and tangents.
    EvaluationData init_evaluation(
        const CurrentState& state,
        bool                with_strain,
        bool                with_B,
        bool                with_G,
        bool                with_resultants,
        const Field*        ip_stress    = nullptr,
        int                 ip_start_idx = 0
    ) const;
    void load_ip_resultants(EvaluationData& data, const Field& ip_stress, int ip_start_idx) const;
    void compute_material_resultants(EvaluationData& data) const;

    // Section and material stiffness helpers.
    // These helpers isolate access to the assigned ShellSection and the
    // topology-density stiffness scale. resultant_stiffness() evaluates the
    // zero-strain generalized tangent at a physical reference point and applies
    // the same scale used by nonlinear section resultants.
    ShellSection* shell_section() const;
    Precision topology_stiffness_scale() const;
    Mat8 resultant_stiffness(Precision r = Precision(0), Precision s = Precision(0)) const;

    // Element tangent and force assembly internals.
    // These routines assemble the pieces of the shell residual linearization:
    // constitutive material stiffness, stress-dependent geometric stiffness,
    // internal force, and objective drilling stabilization. The weighted
    // Hessian helper is kept with this group because it is only used while
    // contracting generalized resultants into the geometric tangent.
    void assemble_material_stiffness(const EvaluationData& data, Mat6N& Kmat) const;
    void assemble_geometric_stiffness(const EvaluationData& data, Mat6N& Kgeo) const;
    void assemble_internal_force(const EvaluationData& data, Vec6N& internal_force) const;
    void assemble_drill_stabilization(
        const EvaluationData& data,
        Mat6N*                stiffness_matrix,
        Vec6N*                internal_force
    ) const;
    void add_weighted_natural_hessian(
        const EvaluationData& data,
        const ReferencePoint& point,
        const Vec8&           weights,
        Mat6N&                Kgeo
    ) const;

    // Generalized and physical result recovery helpers.
    // Output routines need values at arbitrary natural coordinates, not only at
    // integration points. This group reconstructs generalized strain/resultant
    // values, evaluates the shell deformation gradient through the thickness,
    // and converts generalized section response into physical six-component
    // stress and strain vectors.
    Vec8 generalized_strain_at(
        const EvaluationData& data,
        const Vec6N&          q,
        Precision             r,
        Precision             s,
        bool                  nonlinear
    ) const;
    Vec8 generalized_resultant_at(
        const EvaluationData& data,
        const Vec6N&          q,
        Precision             r,
        Precision             s,
        bool                  nonlinear,
        Vec8*                 strain_out = nullptr
    ) const;
    Mat3 deformation_gradient_at(
        const CurrentState& state,
        Precision           r,
        Precision           s,
        Precision           z
    ) const;
    void physical_stress_strain_at(
        const EvaluationData& data,
        const Vec6N&          q,
        Precision             r,
        Precision             s,
        Precision             zeta,
        bool                  nonlinear,
        Vec6&                 strain_out,
        Vec6&                 stress_out
    ) const;

    // Structural stiffness interface.
    // These are the ElementStructural callbacks used by the solvers. stiffness()
    // returns the linear/material shell tangent, stiffness_geom() rebuilds the
    // stress-dependent tangent from stored integration-point resultants, and
    // stiffness_tangent() performs the full nonlinear update including material,
    // geometric and drilling terms plus nodal internal-force scattering.
    MapMatrix stiffness(Precision* buffer) override;
    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override;
    MapMatrix stiffness_tangent(
        Precision*   buffer,
        Field&       ip_stress_state,
        NodeData&    nodal_forces,
        const Field& displacement
    ) override;

    // Basic geometric integration and mass interface.
    // These functions expose scalar geometric properties used outside the
    // nonlinear shell tangent path. They integrate over the curved reference
    // midsurface and assemble consistent translational and rotational mass terms.
    Precision volume() override;
    MapMatrix mass(Precision* buffer) override;

    // Volume field integration interface.
    // These callbacks integrate user-supplied scalar, vector and tensor fields
    // through the shell reference volume. Vector fields can either be returned
    // as one element resultant or consistently distributed to nodal load
    // entries. The fields are evaluated on the current midsurface while the
    // integration measure remains the reference volume measure.
    Precision integrate_scalar_field(bool               scale_by_density,
                                     const ScalarField& field) override;
    Vec3      integrate_vector_field(bool            scale_by_density,
                                     const VecField& field) override;
    void      integrate_vector_field(Field&          node_loads,
                                     bool            scale_by_density,
                                     const VecField& field) override;
    Mat3      integrate_tensor_field(bool            scale_by_density,
                                     const TenField& field) override;

    // Result, force and compatibility callbacks from the structural interface.
    // This group contains post-processing and solver helper functions that have
    // to match the behavior of the other element families. Shell-specific
    // stress states and section forces are recovered where possible; generic
    // beam-style shear flow and beam section force callbacks report that they
    // are not available for shell elements.
    void compute_stress_strain(
        Field*           strain,
        Field*           stress,
        const Field&     displacement,
        const RowMatrix& rst,
        int              offset,
        bool             use_green_lagrange_nl
    ) override;
    void compute_stress_state(
        Field&       stress_state,
        const Field& displacement,
        int          offset,
        bool         use_green_lagrange_nl
    ) override;
    void compute_internal_force_nonlinear(
        Field&       node_forces,
        const Field& ip_stress
    ) override;
    void compute_compliance(
        Field& displacement,
        Field& result
    ) override;
    void compute_compliance_angle_derivative(
        Field& displacement,
        Field& result
    ) override;
    bool compute_shear_flow(
        Field&       shear_flow,
        const Field& displacement,
        int          offset
    ) override;
    bool compute_beam_section_forces(
        Field&       section_forces,
        const Field& displacement,
        int          offset
    ) override;
    bool compute_shell_section_forces(
        Field&       section_forces,
        Field&       contribution_count,
        const Field& displacement
    ) override;

};

// Suppress implicit template generation outside the dedicated shell translation
// unit, which explicitly instantiates all supported topologies.
extern template struct FRTShell<3>;
extern template struct FRTShell<4>;
extern template struct FRTShell<6>;
extern template struct FRTShell<8>;

} // namespace fem::model
