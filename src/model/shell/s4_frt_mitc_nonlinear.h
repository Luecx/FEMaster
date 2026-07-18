    /**
     * @file s4_frt_mitc_nonlinear.h
     * @brief Four-node fully geometrically nonlinear FRT/MITC shell element.
     *
     * The element is a Total-Lagrangian Reissner-Mindlin shell with finite nodal
     * rotations. Each solver node contributes six global displacement coordinates
     * ordered as [ux, uy, uz, rx, ry, rz]. The rotational coordinates form one
     * total global axis-angle vector theta. They are updated additively by the
     * nonlinear solver, while the physical director is evaluated exactly from
     * R(theta) = Exp(skew(theta)). The element tangent is differentiated exactly
     * with respect to the three total rotation-vector components.
     *
     * Mechanical model:
     *   - Total-Lagrangian reference-midsurface integration.
     *   - Current director d = R d0 with R = Exp(theta).
     *   - Green-Lagrange membrane, bending and transverse-shear strains.
     *   - MITC4 transverse-shear tying, including differentiated B and G terms.
     *   - Optional compact membrane EAS condensation controlled by eas_parameters.
     *   - Split nonlinear tangent K_T = K_mat + K_geo.
     *   - K_geo contains membrane, bending and transverse-shear resultant terms.
     *   - compute_stress_state stores the eight generalized shell resultants
     *     [N11,N22,N12,M11,M22,M12,Q13,Q23] at the integration points.
     *   - compute_stress_strain reconstructs physical through-thickness stress and
     *     strain tensors at the requested normalized thickness coordinate t.
     *
     * A small consistent drilling penalty regularizes the otherwise redundant
     * rotation about the reference director. Its force and stiffness are both
     * included, so the assembled Newton tangent remains the derivative of the
     * internal-force vector.
     *
     * The reference director and in-plane basis are generated locally from the
     * planar reference element geometry.
     */

    #pragma once

    #include "shell.h"
    #include "../geometry/surface/surface4.h"

    #include <array>
    #include <memory>
    #include <string>

    namespace fem::model {

    struct S4FRTMITC : ShellElement<4> {
        using Base = ShellElement<4>;

        static constexpr Index num_nodes      = 4;
        static constexpr Index dofs_per_node = 6;
        static constexpr Index num_dofs       = num_nodes * dofs_per_node;
        static constexpr Index num_strains    = 8;

        // 0 disables EAS. Use 4 for the compact membrane EAS4 mode or 5 for EAS5.
        static constexpr Index eas_parameters = 0;

        // Small consistent drill stabilization. The penalty acts on the total
        // rotation-vector component along the reference director. It is kept small
        // and is not part of the eight physical shell resultants.
        static constexpr Precision drill_scale = Precision(1e-5);

        using Vec8     = StaticVector<num_strains>;
        using Vec24    = StaticVector<num_dofs>;
        using Mat8     = StaticMatrix<num_strains, num_strains>;
        using Mat24    = StaticMatrix<num_dofs, num_dofs>;
        using Mat8x24  = StaticMatrix<num_strains, num_dofs>;
        using Mat43    = StaticMatrix<4, 3>;
        using Mat42    = StaticMatrix<4, 2>;
        using Vec24Mat = std::array<Mat24, num_strains>;

        struct VectorDerivatives {
            Vec3                        value = Vec3::Zero();
            StaticMatrix<num_dofs, 3>   d1    = StaticMatrix<num_dofs, 3>::Zero();
            std::array<Mat24, 3>        d2{};

            VectorDerivatives() {
                for (auto& matrix : d2) {
                    matrix.setZero();
                }
            }
        };

        struct ScalarDerivatives {
            Precision value = Precision(0);
            Vec24     d1    = Vec24::Zero();
            Mat24     d2    = Mat24::Zero();
        };

        struct CurrentState {
            Mat43 x     = Mat43::Zero();
            Mat43 d     = Mat43::Zero();
            Mat43 theta = Mat43::Zero();
        };

        struct EvaluationData {
            // -----------------------------------------------------------------
            // Requested evaluation content
            // -----------------------------------------------------------------

            bool with_strain     = false;
            bool with_B          = false;
            bool with_G          = false;
            bool with_resultants = false;
            bool with_eas        = false;

            // -----------------------------------------------------------------
            // Current nodal state and derivative data
            // -----------------------------------------------------------------

            CurrentState                             state;
            std::array<VectorDerivatives, num_nodes> x_nodes;
            std::array<VectorDerivatives, num_nodes> d_nodes;

            // -----------------------------------------------------------------
            // Generalized shell stiffness/resultants
            // -----------------------------------------------------------------

            Mat8 H = Mat8::Zero();

            // -----------------------------------------------------------------
            // MITC4 tying data
            //
            // Index convention:
            //   0 : bottom tying point, r =  0, s = -1
            //   1 : top    tying point, r =  0, s =  1
            //   2 : left   tying point, r = -1, s =  0
            //   3 : right  tying point, r =  1, s =  0
            // -----------------------------------------------------------------

            std::array<Vec2,                      4> tying_shear_nat;
            std::array<StaticMatrix<2, num_dofs>, 4> tying_B_nat;
            std::array<std::array<Mat24, 2>,      4> tying_G_nat;

            // -----------------------------------------------------------------
            // Integration point data
            //
            // Index convention:
            //   0..3 : Gauss integration points.
            //
            // ip_strain:
            //   [eps11, eps22, gamma12, kappa11, kappa22, kappa12,
            //    gamma13, gamma23]
            //
            // ip_resultants:
            //   [N11, N22, N12, M11, M22, M12, Q13, Q23]
            // -----------------------------------------------------------------

            std::array<Vec8,      4> ip_strain;
            std::array<Mat8x24,   4> ip_B;
            std::array<Vec24Mat,  4> ip_G;
            std::array<Vec8,      4> ip_resultants;
            std::array<Mat8,      4> ip_tangent;
            std::array<Precision, 4> ip_weight;

            // -----------------------------------------------------------------
            // EAS condensation data
            // -----------------------------------------------------------------

            StaticMatrix<eas_parameters, eas_parameters> eas_Kaa;
            StaticMatrix<eas_parameters, num_dofs>       eas_Jb;
            StaticVector<eas_parameters>                 eas_ba;
            StaticVector<eas_parameters>                 eas_alpha;

            // -----------------------------------------------------------------
            // Drill stabilization
            // -----------------------------------------------------------------

            Precision drill_k = Precision(0);

            EvaluationData();
        };

        struct ReferenceData {
            // -----------------------------------------------------------------
            // Reference nodal geometry
            // -----------------------------------------------------------------

            Mat43 X  = Mat43::Zero();
            Mat43 d0 = Mat43::Zero();
            Mat42 xy = Mat42::Zero();

            Vec3 e1 = Vec3::Zero();
            Vec3 e2 = Vec3::Zero();
            Vec3 e3 = Vec3::Zero();

            Precision area = Precision(0);

            // -----------------------------------------------------------------
            // Reference point indexing
            //
            //   0..3 : Gauss integration points
            //   4    : bottom tying point, r =  0, s = -1
            //   5    : top    tying point, r =  0, s =  1
            //   6    : left   tying point, r = -1, s =  0
            //   7    : right  tying point, r =  1, s =  0
            //
            // A temporary ReferenceData object may also be used with ref_id = 0
            // for arbitrary output points that are not part of the cached eight
            // points above.
            // -----------------------------------------------------------------

            static constexpr Index num_ip                 = 4;
            static constexpr Index num_tying_points       = 4;
            static constexpr Index tying_start            = num_ip;
            static constexpr Index num_ref_points         = num_ip
                                                           + num_tying_points;

            StaticVector<num_ref_points> r    = StaticVector<num_ref_points>::Zero();
            StaticVector<num_ref_points> s    = StaticVector<num_ref_points>::Zero();
            StaticVector<num_ref_points> w    = StaticVector<num_ref_points>::Zero();
            StaticVector<num_ref_points> detJ = StaticVector<num_ref_points>::Zero();

            std::array<StaticVector<4>,    num_ref_points> N{};
            std::array<StaticMatrix<4, 2>, num_ref_points> dN_rs{};

            std::array<StaticVector<4>, num_ref_points> dN_da{};
            std::array<StaticVector<4>, num_ref_points> dN_db{};

            std::array<Mat2, num_ref_points> A{};
            std::array<Mat2, num_ref_points> invA{};

            // -----------------------------------------------------------------
            // Reference midsurface/director derivatives at all reference points
            // -----------------------------------------------------------------

            std::array<Vec3, num_ref_points> X_a{};
            std::array<Vec3, num_ref_points> X_b{};
            std::array<Vec3, num_ref_points> X_xi{};
            std::array<Vec3, num_ref_points> X_eta{};

            std::array<Vec3, num_ref_points> D{};
            std::array<Vec3, num_ref_points> D_a{};
            std::array<Vec3, num_ref_points> D_b{};

            ReferenceData();
        };

        Surface4                        geometry;
        math::quadrature::Quadrature          integration_scheme_;
        std::unique_ptr<ReferenceData>  reference_data_;

        S4FRTMITC(ID id, std::array<ID, 4> nodes);

        ~S4FRTMITC() override = default;

        void step_begin() override;
        void step_end  () override;

        Mat43              init_ref_node_coords() const;
        void init_reference_basis(ReferenceData& data) const;
        void init_ref_point_data(ReferenceData& data, Index ref_id, Precision r, Precision s, Precision w) const;
        std::string type_name() const override;

        std::shared_ptr<SurfaceInterface> surface(int surface_id) override;

        const math::quadrature::Quadrature& integration_scheme() const override;

        StaticVector<4> shape_function(Precision r, Precision s) const;

        StaticMatrix<4, 2> shape_derivative(Precision r, Precision s) const;

        ShellSection* shell_section() const;

        // ---------------------------------------------------------------------
        // Basic vector/tensor helpers
        // ---------------------------------------------------------------------

        static ScalarDerivatives dot_derivatives(const VectorDerivatives& a,
                                                 const VectorDerivatives& b);

        static ScalarDerivatives scaled(const ScalarDerivatives& a, Precision scale);

        static VectorDerivatives linear_combination(const std::array<VectorDerivatives, 4>& values,
                                                    const StaticVector<4>&                  coeffs);

        // ---------------------------------------------------------------------
        // Reference/current kinematics
        // ---------------------------------------------------------------------


        StaticMatrix<4, 6> node_coords_current_6() const;

        const ReferenceData& reference_data() const;

        int cached_reference_point_id(Precision r, Precision s) const;

        CurrentState current_state() const;

        CurrentState reference_state() const;

        CurrentState current_state_from_displacement(
            const Field&        displacement) const;

        Vec24 element_displacement_vector(const Field& displacement) const;

        Mat3 reference_basis_global() const;

        void shape_gradients_physical(const StaticMatrix<4, 2>& dN_rs,
                                      StaticVector<4>&          dN_da,
                                      StaticVector<4>&          dN_db,
                                      Precision&                detJ,
                                      Mat2&                     A) const;

        std::array<VectorDerivatives, 4> x_derivatives(const CurrentState& state) const;

        std::array<VectorDerivatives, 4> director_derivatives(
            const CurrentState& state,
            bool                with_second_derivatives = true) const;

        EvaluationData init_evaluation(
            const CurrentState& state,
            bool                with_strain,
            bool                with_B,
            bool                with_G,
            bool                with_resultants,
            bool                with_eas,
            const Field*        ip_stress    = nullptr,
            int                 ip_start_idx = 0
        ) const;

        void reference_fields(const StaticVector<4>& N,
                              const StaticVector<4>& dN_da,
                              const StaticVector<4>& dN_db,
                              Vec3& X_a,
                              Vec3& X_b,
                              Vec3& D,
                              Vec3& D_a,
                              Vec3& D_b) const;

        void compute_raw_strain(
            const EvaluationData& data,
            const ReferenceData&  ref,
            Index                 ref_id,
            Vec8&                 strain,
            Mat8x24*                  B = nullptr,
            Vec24Mat*                 G = nullptr
        ) const;

        void compute_natural_shear(
            const EvaluationData&      data,
            const ReferenceData&       ref,
            Index                      ref_id,
            Vec2&                      shear_nat,
            StaticMatrix<2, num_dofs>*       B_nat = nullptr,
            std::array<Mat24, 2>*            G_nat = nullptr
        ) const;

        void compute_mitc4_shear(
            const EvaluationData&      data,
            const ReferenceData&       ref,
            Index                      ref_id,
            Vec2&                      shear,
            StaticMatrix<2, num_dofs>*       B_shear = nullptr,
            std::array<Mat24, 2>*            G_shear = nullptr
        ) const;

        void evaluate_tying_point(
            EvaluationData&      data,
            Index                tying_id,
            const ReferenceData& ref,
            Index                ref_id
        ) const;

        void evaluate_integration_point(
            EvaluationData&      data,
            Index                ip,
            const ReferenceData& ref,
            Index                ref_id
        ) const;

        void load_ip_resultants(
            EvaluationData& data,
            const Field&    ip_stress,
            int             ip_start_idx
        ) const;

        void compute_eas_data(EvaluationData& data) const;
        void compute_material_resultants(EvaluationData& data) const;

        void assemble_material_stiffness(
            const EvaluationData& data,
            Mat24&                Kmat
        ) const;

        void assemble_geometric_stiffness(
            const EvaluationData& data,
            Mat24&                Kgeo
        ) const;

        void assemble_internal_force(
            const EvaluationData& data,
            Vec24&                internal_force
        ) const;

        // ---------------------------------------------------------------------
        // Material/resultant matrices and EAS
        // ---------------------------------------------------------------------

        Precision topology_stiffness_scale() const;
        Mat8 resultant_stiffness() const;

        StaticMatrix<num_strains, eas_parameters> eas_matrix(Precision r, Precision s) const;

        StaticVector<eas_parameters> compute_eas_alpha(
            const CurrentState& state,
            const Mat8&         H
        ) const;

        StaticVector<eas_parameters> compute_eas_alpha_linear(
            const EvaluationData& data,
            const Vec24&          q
        ) const;

        Vec8 generalized_strain_at(
            const EvaluationData&                 data,
            const Vec24&                          q,
            const StaticVector<eas_parameters>&   alpha,
            Precision                             r,
            Precision                             s,
            bool                                  nonlinear
        ) const;

        Vec8 generalized_resultant_at(
            const EvaluationData&                 data,
            const Vec24&                          q,
            const StaticVector<eas_parameters>&   alpha,
            Precision                             r,
            Precision                             s,
            bool                                  nonlinear,
            Vec8*                                 strain_out = nullptr
        ) const;

        Mat3 current_output_basis(const CurrentState& state,
                                  Precision           r,
                                  Precision           s) const;

        Mat3 deformation_gradient_at(const CurrentState& state,
                                     Precision           r,
                                     Precision           s,
                                     Precision           z) const;

        void physical_stress_strain_at(
            const EvaluationData&                 data,
            const Vec24&                          q,
            const StaticVector<eas_parameters>&   alpha,
            Precision                             r,
            Precision                             s,
            Precision                             zeta,
            bool                                  nonlinear,
            Vec6&                                 strain_out,
            Vec6&                                 stress_out
        ) const;

        // ---------------------------------------------------------------------
        // Drill stabilization and StructuralElement interface
        // ---------------------------------------------------------------------

        Precision drill_stiffness_per_node(const Mat8& H) const;

        void add_drill_stiffness(Mat24&              stiffness_matrix,
                                 const Mat8&         H) const;

        void add_drill_force(Vec24&              internal_force,
                             const CurrentState& state,
                             const Mat8&         H) const;

        MapMatrix stiffness(Precision* buffer) override;

        MapMatrix stiffness_geom(Precision*   buffer,
                                 const Field& ip_stress,
                                 int          ip_start_idx) override;

        MapMatrix stiffness_tangent(Precision* buffer,
                                    Field&       ip_stress_state,
                                    NodeData&    nodal_forces,
                                    const Field& displacement) override;

        RowMatrix stress_strain_ip_rst() override;

        RowMatrix stress_strain_nodal_rst() override;

        // ---------------------------------------------------------------------
        // Mass, volume and load integration helpers
        // ---------------------------------------------------------------------

        Precision volume() override;

        StaticMatrix<4, 4> integrate_NNt() const;

        MapMatrix mass(Precision* buffer) override;

        Vec3 global_point_current(Precision r, Precision s) const;

        // integration over the volume
        Precision integrate_scalar_field(bool            scale_by_density,
                                         const ScalarField& field) override;
        Vec3      integrate_vector_field(bool            scale_by_density,
                                         const VecField& field) override;
        void      integrate_vector_field(Field&          node_loads,
                                         bool            scale_by_density,
                                         const VecField& field) override;
        Mat3      integrate_tensor_field(bool            scale_by_density,
                                         const TenField& field) override;

        // compute functions (inherited)
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

    } // namespace fem::model
