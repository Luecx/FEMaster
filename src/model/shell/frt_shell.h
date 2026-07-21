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
 * The first strain derivative is denoted by B and the eight second strain
 * derivatives are denoted by G:
 *
 *     B    = d(epsilon) / d(q),
 *     G[a] = d²(epsilon[a]) / d(q)².
 *
 * These quantities yield the consistent internal force and tangent
 *
 *     f_int = integral_A0 B^T n dA0,
 *     K_mat = integral_A0 B^T H B dA0,
 *     K_geo = integral_A0 sum_a n[a] G[a] dA0.
 *
 * EAS and follower-load contributions intentionally do not belong to this
 * implementation. The concrete MITC elements only replace the selected
 * compatible strain components and apply exactly the same linear interpolation
 * to their values, B rows and G matrices.
 *
 * @see FRTShellS3
 * @see FRTShellS4
 * @see FRTShellS6
 * @see FRTShellS8
 *
 * @author Finn Eggers
 * @date 20.07.2026
 */

#pragma once

#include "shell.h"

#include "../../math/quadrature.h"

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace fem::model {

/**
 * @brief Common finite-rotation Total-Lagrangian Reissner-Mindlin shell.
 *
 * `FRTShell<N>` contains the nonlinear nodal director kinematics, the exact
 * first and second strain derivatives, the pointwise reference-surface
 * geometry, constitutive section evaluation, element assembly, drilling
 * stabilization, mass integration and result recovery shared by all supported
 * shell topologies.
 *
 * Every node carries the global degrees of freedom
 * `[ux, uy, uz, rx, ry, rz]`. The rotational coordinates form one total
 * axis-angle vector. The current nodal director is evaluated from
 * `d = Exp(theta) d0`, including exact first and second derivatives with
 * respect to the three rotational coordinates.
 *
 * Concrete elements must provide:
 *
 * - shape functions and first natural derivatives,
 * - natural nodal coordinates,
 * - numerical quadrature and surface connectivity,
 * - MITC tying-point coordinates,
 * - interpolation of the selected assumed natural strain components.
 *
 * The common class stores the actual dimensions in its matrix aliases. The
 * number of integration and tying points remains element-specific and is held
 * in vectors initialized once at `step_begin()`.
 *
 * @tparam N Number of shell midsurface nodes.
 */
template<Index N>
struct FRTShell : ShellElement<N> {
    using Base = ShellElement<N>;

    static_assert(N == 3 || N == 4 || N == 6 || N == 8,
                  "FRTShell supports only S3, S4, S6 and S8 topologies");

    static constexpr Index num_nodes     = N;
    static constexpr Index dofs_per_node = 6;
    static constexpr Index num_dofs      = num_nodes * dofs_per_node;
    static constexpr Index num_strains   = 8;

    // Small consistent drilling penalty. The penalty acts on the component of
    // the total rotation vector along the nodal reference director.
    static constexpr Precision drill_scale = Precision(1e-5);

    // Fixed-size element vectors and matrices
    using Vec8      = StaticVector<num_strains>;
    using Vec6N     = StaticVector<num_dofs>;
    using Mat8      = StaticMatrix<num_strains, num_strains>;
    using Mat6N     = StaticMatrix<num_dofs, num_dofs>;
    using Mat8x6N   = StaticMatrix<num_strains, num_dofs>;
    using Mat2x6N   = StaticMatrix<2, num_dofs>;
    using Mat3x6N   = StaticMatrix<3, num_dofs>;
    using MatN3     = StaticMatrix<num_nodes, 3>;
    using MatN2     = StaticMatrix<num_nodes, 2>;
    using MatN6     = StaticMatrix<num_nodes, 6>;
    using VecN      = StaticVector<num_nodes>;
    using Vec6NMat  = std::array<Mat6N, num_strains>;
    using Vec3G     = std::array<Mat6N, 3>;

    /**
     * @brief Value and derivatives of a three-dimensional vector quantity.
     *
     * The derivative convention is
     *
     *     d1(i,c)    = d(value[c]) / d(q[i]),
     *     d2[c](i,j) = d²(value[c]) / (d(q[i]) d(q[j])).
     *
     * It is used for current positions, rotated directors and their spatial
     * derivatives. The second derivatives are required by the consistent
     * geometric tangent.
     */
    struct VectorDerivatives {
        Vec3                      value = Vec3::Zero();
        StaticMatrix<num_dofs, 3> d1    = StaticMatrix<num_dofs, 3>::Zero();
        Vec3G                     d2{};

        VectorDerivatives();
    };

    /**
     * @brief Value, gradient and Hessian of a scalar kinematic quantity.
     *
     * Scalar products such as `x_,a . x_,b`, `x_,a . d_,b` and
     * `x_,a . d` are represented by this type while constructing membrane,
     * bending and transverse-shear strains.
     */
    struct ScalarDerivatives {
        Precision value = Precision(0);
        Vec6N     d1    = Vec6N::Zero();
        Mat6N     d2    = Mat6N::Zero();
    };

    /**
     * @brief Current nodal positions, directors and total rotation vectors.
     *
     * `x` contains the current midsurface node positions, `d` contains the
     * rotated nodal directors and `theta` contains the global total axis-angle
     * rotation vectors.
     */
    struct CurrentState {
        MatN3 x     = MatN3::Zero();
        MatN3 d     = MatN3::Zero();
        MatN3 theta = MatN3::Zero();
    };

    /**
     * @brief Reference geometry evaluated at one natural shell point.
     *
     * The point stores the shape functions, the actual curved reference
     * tangents, the complete surface Jacobian and the local orthonormal basis.
     * `A` maps local tangent derivatives `(a,b)` to natural derivatives `(r,s)`:
     *
     *     [N_,r, N_,s]^T = A [N_,a, N_,b]^T.
     *
     * The reference director field and its derivatives are interpolated from
     * the nodal reference directors using the same shell shape functions.
     */
    struct ReferencePoint {
        Precision r      = Precision(0);
        Precision s      = Precision(0);
        Precision w      = Precision(0);
        Precision detJ   = Precision(0);

        VecN  shape      = VecN::Zero();
        MatN2 dshape_rs  = MatN2::Zero();
        VecN  dshape_a   = VecN::Zero();
        VecN  dshape_b   = VecN::Zero();

        Vec3 X_r = Vec3::Zero();
        Vec3 X_s = Vec3::Zero();
        Vec3 X_a = Vec3::Zero();
        Vec3 X_b = Vec3::Zero();

        Vec3 D   = Vec3::Zero();
        Vec3 D_r = Vec3::Zero();
        Vec3 D_s = Vec3::Zero();
        Vec3 D_a = Vec3::Zero();
        Vec3 D_b = Vec3::Zero();

        Vec3 e1 = Vec3::Zero();
        Vec3 e2 = Vec3::Zero();
        Vec3 e3 = Vec3::Zero();

        Mat2 A    = Mat2::Zero();
        Mat2 invA = Mat2::Zero();
    };

    /**
     * @brief Persistent reference configuration of one shell element.
     *
     * The data is initialized at `step_begin()` and remains unchanged during
     * all Newton iterations of the step. Integration and tying points both use
     * the actual isoparametric reference surface.
     */
    struct ReferenceData {
        // Reference nodal midsurface positions and directors
        MatN3 X  = MatN3::Zero();
        MatN3 d0 = MatN3::Zero();

        // Cached integration and MITC tying points
        std::vector<ReferencePoint> ip_points;
        std::vector<ReferencePoint> tying_points;

        // Reference midsurface area obtained from the integration rule
        Precision area = Precision(0);
    };

    /**
     * @brief Complete kinematic and constitutive data for one element evaluation.
     *
     * The compatible natural strains are evaluated at every tying point. The
     * concrete element then interpolates its assumed strain field at each
     * integration point. The resulting local generalized strains, B matrices,
     * G matrices, section resultants and tangents are stored per integration
     * point.
     */
    struct EvaluationData {
        // Requested evaluation content. Computing G implies B and strain;
        // computing B implies strain.
        bool with_strain     = false;
        bool with_B          = false;
        bool with_G          = false;
        bool with_resultants = false;

        // Current nodal state and exact position/director derivatives
        CurrentState                             state;
        std::array<VectorDerivatives, num_nodes> x_nodes;
        std::array<VectorDerivatives, num_nodes> d_nodes;

        // Homogeneous section tangent used for drilling stabilization
        Mat8 H = Mat8::Zero();

        // Compatible generalized strains at MITC tying points in natural
        // coordinates [r,s]. Engineering shear components are used.
        std::vector<Vec8>     tying_strain_nat;
        std::vector<Mat8x6N>  tying_B_nat;
        std::vector<Vec6NMat> tying_G_nat;

        // Generalized integration-point kinematics in the local orthonormal
        // reference basis. ip_B = d(ip_strain)/d(q) and
        // ip_G[a] = d²(ip_strain[a])/d(q)².
        std::vector<Vec8>     ip_strain;
        std::vector<Mat8x6N>  ip_B;
        std::vector<Vec6NMat> ip_G;

        // Integration-point section response. ip_tangent is the local section
        // tangent H = d(ip_resultants)/d(ip_strain), not the element tangent.
        std::vector<Vec8> ip_resultants;
        std::vector<Mat8> ip_tangent;

        // Complete reference-area quadrature weights w * detJ
        std::vector<Precision> ip_weight;

        // Consistent drilling stabilization stiffness per node
        Precision drill_k = Precision(0);

        EvaluationData(Index num_ip, Index num_tying);
    };

    // Persistent reference state initialized for the current analysis step
    std::unique_ptr<ReferenceData> reference_data_;

    // Construction and step state
    FRTShell(ID id, const std::array<ID, N>& nodes);
    ~FRTShell() override = default;

    void step_begin() override;
    void step_end  () override;

    // Concrete interpolation and MITC formulation. Natural coordinates and
    // node ordering must agree with the associated Surface implementation.
    virtual VecN  shape_function         (Precision r, Precision s) const = 0;
    virtual MatN2 shape_derivative       (Precision r, Precision s) const = 0;
    virtual MatN2 node_coords_natural    () const = 0;
    virtual std::vector<Vec2> tying_point_coordinates() const = 0;

    // The concrete MITC implementation receives compatible generalized strains
    // in natural coordinates and replaces only the components belonging to its
    // assumed-strain space. Values, B rows and G matrices must be interpolated
    // with identical coefficients.
    virtual void apply_mitc_natural(
        const EvaluationData& data,
        const ReferencePoint& point,
        Vec8&                 strain_nat,
        Mat8x6N*              B_nat,
        Vec6NMat*             G_nat
    ) const = 0;

    // Reference configuration and pointwise curved-surface geometry
    const ReferenceData& reference_data() const;
    MatN3                init_reference_node_coords() const;
    MatN3                init_reference_directors(const MatN3& X) const;
    ReferencePoint       make_reference_point(Precision r, Precision s, Precision w) const;
    Vec3                 reference_position(Precision r, Precision s) const;
    Mat3                 reference_basis_global(Precision r, Precision s) const;
    const ReferencePoint* cached_reference_point(Precision r, Precision s) const;

    // Current nodal configuration and exact director derivatives
    MatN6        node_coords_current_6() const;
    CurrentState current_state() const;
    CurrentState reference_state() const;
    CurrentState current_state_from_displacement(const Field& displacement) const;
    Vec6N        element_displacement_vector(const Field& displacement) const;

    std::array<VectorDerivatives, N> x_derivatives(const CurrentState& state) const;
    std::array<VectorDerivatives, N> director_derivatives(
        const CurrentState& state,
        bool                with_second_derivatives = true
    ) const;

    // Scalar and vector derivative operations used by the shell strains
    static ScalarDerivatives dot_derivatives(const VectorDerivatives& a,
                                             const VectorDerivatives& b);
    static ScalarDerivatives scaled(const ScalarDerivatives& value, Precision scale);
    static VectorDerivatives linear_combination(
        const std::array<VectorDerivatives, N>& values,
        const VecN&                            coefficients
    );

    // Compatible Total-Lagrangian generalized shell strains in natural
    // coordinates. The ordering matches the local generalized strain vector,
    // but components refer to r and s before the pointwise basis transformation.
    void compute_natural_strain(
        const EvaluationData& data,
        const ReferencePoint& point,
        Vec8&                 strain_nat,
        Mat8x6N*              B_nat = nullptr,
        Vec6NMat*             G_nat = nullptr
    ) const;

    // Linear reference-geometry transformations from natural covariant strain
    // components to the local orthonormal tangent basis. The same maps are
    // applied to values, B rows and G matrices.
    void transform_strain_to_local(
        const ReferencePoint& point,
        Vec8&                 strain,
        Mat8x6N*              B = nullptr,
        Vec6NMat*             G = nullptr
    ) const;

    StaticMatrix<3, 3> in_plane_transform(const ReferencePoint& from,
                                          const ReferencePoint& to) const;

    // Complete element evaluation. Tying values are prepared first, followed by
    // MITC interpolation, local-basis transformation and optional section
    // response at every integration point.
    EvaluationData init_evaluation(
        const CurrentState& state,
        bool                with_strain,
        bool                with_B,
        bool                with_G,
        bool                with_resultants,
        const Field*        ip_stress    = nullptr,
        int                 ip_start_idx = 0
    ) const;

    void evaluate_tying_points(EvaluationData& data) const;
    void evaluate_integration_points(EvaluationData& data) const;
    void load_ip_resultants(EvaluationData& data, const Field& ip_stress, int ip_start_idx) const;
    void compute_material_resultants(EvaluationData& data) const;

    // Section response and consistent element assembly
    ShellSection* shell_section() const;
    Precision     topology_stiffness_scale() const;
    Mat8          resultant_stiffness(Precision r = Precision(0), Precision s = Precision(0)) const;

    void assemble_material_stiffness(const EvaluationData& data, Mat6N& Kmat) const;
    void assemble_geometric_stiffness(const EvaluationData& data, Mat6N& Kgeo) const;
    void assemble_internal_force    (const EvaluationData& data, Vec6N& internal_force) const;

    // Generalized and physical result recovery at arbitrary natural points
    Vec8 generalized_strain_at(const EvaluationData& data,
                               const Vec6N&           q,
                               Precision              r,
                               Precision              s,
                               bool                   nonlinear) const;

    Vec8 generalized_resultant_at(const EvaluationData& data,
                                  const Vec6N&           q,
                                  Precision              r,
                                  Precision              s,
                                  bool                   nonlinear,
                                  Vec8*                  strain_out = nullptr) const;

    Mat3 current_output_basis(const CurrentState& state, Precision r, Precision s) const;
    Mat3 deformation_gradient_at(const CurrentState& state, Precision r, Precision s, Precision z) const;

    void physical_stress_strain_at(const EvaluationData& data,
                                   const Vec6N&           q,
                                   Precision              r,
                                   Precision              s,
                                   Precision              zeta,
                                   bool                   nonlinear,
                                   Vec6&                  strain_out,
                                   Vec6&                  stress_out) const;

    // Consistent drilling stabilization. Rotation about the nodal reference
    // director is regularized without contributing to the physical shell
    // resultants.
    Precision drill_stiffness_per_node(const Mat8& H) const;
    void      add_drill_stiffness(Mat6N& stiffness_matrix, const Mat8& H) const;
    void      add_drill_force(Vec6N& internal_force, const CurrentState& state, const Mat8& H) const;

    // Structural element stiffness and nonlinear internal force
    MapMatrix stiffness(Precision* buffer) override;
    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override;
    MapMatrix stiffness_tangent(Precision* buffer,
                                Field&     ip_stress_state,
                                NodeData&  nodal_forces,
                                const Field& displacement) override;

    void compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) override;

    // Stress, strain and generalized shell-resultant output
    void compute_stress_strain(Field*           strain,
                               Field*           stress,
                               const Field&     displacement,
                               const RowMatrix& rst,
                               int              offset,
                               bool             use_green_lagrange_nl) override;

    void compute_stress_state(Field&       stress_state,
                              const Field& displacement,
                              int          offset,
                              bool         use_green_lagrange_nl) override;

    bool compute_shell_section_forces(Field&       section_forces,
                                      Field&       contribution_count,
                                      const Field& displacement) override;

    // Reference volume, consistent mass and current field evaluation
    Precision volume() override;
    StaticMatrix<N, N> integrate_NNt() const;
    MapMatrix mass(Precision* buffer) override;
    Vec3 global_point_current(Precision r, Precision s) const;

    // Integration over the shell volume. These loads are body-attached or
    // spatially prescribed fields evaluated on the current midsurface, while
    // the integration measure remains the reference volume measure.
    Precision integrate_scalar_field(bool               scale_by_density,
                                     const ScalarField& field) override;
    Vec3      integrate_vector_field(bool            scale_by_density,
                                     const VecField& field) override;
    void      integrate_vector_field(Field&          node_loads,
                                     bool            scale_by_density,
                                     const VecField& field) override;
    Mat3      integrate_tensor_field(bool            scale_by_density,
                                     const TenField& field) override;

    // Generic result helpers inherited by the structural element interface
    void compute_compliance(Field& displacement, Field& result) override;
    void compute_compliance_angle_derivative(Field& displacement, Field& result) override;
    bool compute_shear_flow(Field& shear_flow, const Field& displacement, int offset) override;
    bool compute_beam_section_forces(Field& section_forces,
                                     const Field& displacement,
                                     int          offset) override;
};

// The template implementation is explicitly instantiated for the four shell
// topologies in the corresponding implementation files.
extern template struct FRTShell<3>;
extern template struct FRTShell<4>;
extern template struct FRTShell<6>;
extern template struct FRTShell<8>;

} // namespace fem::model
