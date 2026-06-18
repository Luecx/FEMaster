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

    struct ElementBasis {
        Vec3  e1 = Vec3::Zero();
        Vec3  e2 = Vec3::Zero();
        Vec3  e3 = Vec3::Zero();
        Mat43 d0 = Mat43::Zero();
        Mat42 xy = Mat42::Zero();
    };

    struct CurrentState {
        Mat43 x     = Mat43::Zero();
        Mat43 d     = Mat43::Zero();
        Mat43 theta = Mat43::Zero();
    };

    struct RotationCoefficients {
        Precision a     = Precision(0);
        Precision b     = Precision(0);
        Precision a_s   = Precision(0);
        Precision b_s   = Precision(0);
        Precision a_ss  = Precision(0);
        Precision b_ss  = Precision(0);
    };

    struct StrainData {
        Vec8     strain = Vec8::Zero();
        Mat8x24  B      = Mat8x24::Zero();
        Vec24Mat G{};
        Precision detJ  = Precision(0);

        StrainData() {
            for (auto& matrix : G) {
                matrix.setZero();
            }
        }
    };

    Surface4                 geometry;
    quadrature::Quadrature   integration_scheme_;

    S4FRTMITC(ID id, std::array<ID, 4> nodes);

    ~S4FRTMITC() override = default;

    std::string type_name() const override;

    std::shared_ptr<SurfaceInterface> surface(int surface_id) override;

    const quadrature::Quadrature& integration_scheme() const override;

    StaticVector<4> shape_function(Precision r, Precision s) const;

    StaticMatrix<4, 2> shape_derivative(Precision r, Precision s) const;

    ShellSection* shell_section() const;

    // ---------------------------------------------------------------------
    // Basic vector/tensor helpers
    // ---------------------------------------------------------------------

    static Vec3 normalized(Vec3 v, const std::string& name);

    static Mat3 skew(const Vec3& v);

    static RotationCoefficients rotation_coefficients(Precision angle_squared);

    static void rotation_exp_derivatives(
        const Vec3&                                rotation_vector,
        Mat3&                                      rotation,
        std::array<Mat3, 3>&                       first,
        std::array<std::array<Mat3, 3>, 3>&        second);

    static Mat3 rotation_exp(const Vec3& rotation_vector);

    static ScalarDerivatives dot_derivatives(const VectorDerivatives& a,
                                             const VectorDerivatives& b);

    static ScalarDerivatives scaled(const ScalarDerivatives& a, Precision scale);

    static VectorDerivatives linear_combination(const std::array<VectorDerivatives, 4>& values,
                                                const StaticVector<4>&                  coeffs);

    // ---------------------------------------------------------------------
    // Reference/current kinematics
    // ---------------------------------------------------------------------

    Mat43 node_coords_reference_xyz() const;

    StaticMatrix<4, 6> node_coords_current_6() const;

    ElementBasis reference_basis() const;

    static void compute_local_reference_coordinates(const Mat43& X, ElementBasis& basis);

    CurrentState current_state(const ElementBasis& basis) const;

    CurrentState reference_state(const ElementBasis& basis) const;

    CurrentState current_state_from_displacement(
        const ElementBasis& basis,
        const Field&        displacement) const;

    Vec24 element_displacement_vector(const Field& displacement) const;

    static Mat3 reference_basis_global(const ElementBasis& basis);

    void shape_gradients_physical(const StaticMatrix<4, 2>& dN_rs,
                                  const ElementBasis&       basis,
                                  StaticVector<4>&          dN_da,
                                  StaticVector<4>&          dN_db,
                                  Precision&                detJ,
                                  Mat2&                     A) const;

    std::array<VectorDerivatives, 4> x_derivatives(const CurrentState& state) const;

    std::array<VectorDerivatives, 4> director_derivatives(
        const CurrentState& state,
        const ElementBasis& basis) const;

    void reference_fields(const ElementBasis& basis,
                          const StaticVector<4>& N,
                          const StaticVector<4>& dN_da,
                          const StaticVector<4>& dN_db,
                          Vec3& X_a,
                          Vec3& X_b,
                          Vec3& D,
                          Vec3& D_a,
                          Vec3& D_b) const;

    StrainData raw_strain_B_G_at(const CurrentState& state,
                                 const ElementBasis& basis,
                                 Precision           r,
                                 Precision           s) const;

    void raw_natural_shear_B_G_at(const CurrentState& state,
                                  const ElementBasis& basis,
                                  Precision           r,
                                  Precision           s,
                                  Vec2&               shear_nat,
                                  StaticMatrix<2, num_dofs>& B_nat,
                                  std::array<Mat24, 2>& G_nat) const;

    void mitc4_shear_B_G_at(const CurrentState& state,
                            const ElementBasis& basis,
                            Precision           r,
                            Precision           s,
                            const Mat2&         A,
                            Vec2&               shear,
                            StaticMatrix<2, num_dofs>& B_shear,
                            std::array<Mat24, 2>& G_shear) const;

    StrainData strain_B_G_at(const CurrentState& state,
                             const ElementBasis& basis,
                             Precision           r,
                             Precision           s) const;

    // ---------------------------------------------------------------------
    // Material/resultant matrices and EAS
    // ---------------------------------------------------------------------

    Mat8 resultant_stiffness() const;

    StaticMatrix<num_strains, eas_parameters> eas_matrix(Precision r, Precision s) const;

    StaticVector<eas_parameters> compute_eas_alpha(const CurrentState& state,
                                                   const ElementBasis& basis,
                                                   const Mat8&         H) const;

    void material_and_geometric_stiffness(const CurrentState& state,
                                            const ElementBasis& basis,
                                            const Mat8&         H,
                                            Mat24&              Kmat,
                                            Mat24&              Kgeo,
                                            Vec24*              force = nullptr) const;

    Vec8 shell_resultant_at(const CurrentState& state,
                            const ElementBasis& basis,
                            const Mat8&         H,
                            Precision           r,
                            Precision           s,
                            Vec8*               strain_out = nullptr) const;

    StaticVector<eas_parameters> compute_eas_alpha_linear(
        const Field&        displacement,
        const ElementBasis& basis,
        const Mat8&         H) const;

    Vec8 generalized_strain_at(const Field&        displacement,
                               const ElementBasis& basis,
                               Precision           r,
                               Precision           s,
                               bool                nonlinear) const;

    Vec8 generalized_resultant_at(const Field&        displacement,
                                  const ElementBasis& basis,
                                  Precision           r,
                                  Precision           s,
                                  bool                nonlinear,
                                  Vec8*               strain_out = nullptr) const;

    Mat3 current_output_basis(const CurrentState& state,
                              const ElementBasis& basis,
                              Precision           r,
                              Precision           s) const;

    Mat3 deformation_gradient_at(const CurrentState& state,
                                 const ElementBasis& basis,
                                 Precision           r,
                                 Precision           s,
                                 Precision           z) const;

    void physical_stress_strain_at(const Field&        displacement,
                                   const ElementBasis& basis,
                                   Precision           r,
                                   Precision           s,
                                   Precision           zeta,
                                   bool                nonlinear,
                                   Vec6&               strain_out,
                                   Vec6&               stress_out) const;

    // ---------------------------------------------------------------------
    // Drill stabilization and StructuralElement interface
    // ---------------------------------------------------------------------

    Precision reference_area(const ElementBasis& basis) const;

    Precision drill_stiffness_per_node(const ElementBasis& basis,
                                       const Mat8&         H) const;

    void add_drill_stiffness(Mat24&              stiffness_matrix,
                             const ElementBasis& basis,
                             const Mat8&         H) const;

    void add_drill_force(Vec24&              internal_force,
                         const CurrentState& state,
                         const ElementBasis& basis,
                         const Mat8&         H) const;

    MapMatrix stiffness(Precision* buffer) override;

    MapMatrix stiffness_geom(Precision*   buffer,
                             const Field& ip_stress,
                             int          ip_start_idx) override;

    void compute_internal_force_nonlinear(Field&       node_forces,
                                          const Field& ip_stress,
                                          int          ip_offset) override;

    void compute_stress_strain(Field*           strain,
                               Field*           stress_out,
                               const Field&     displacement,
                               const RowMatrix& rst,
                               int              offset,
                               bool             use_green_lagrange_nl) override;

    void compute_stress_state(Field&       stress_state,
                              const Field& displacement,
                              int          offset,
                              bool         use_green_lagrange_nl) override;

    RowMatrix stress_strain_ip_rst() override;

    RowMatrix stress_strain_nodal_rst() override;

    // ---------------------------------------------------------------------
    // Mass, volume and load integration helpers
    // ---------------------------------------------------------------------

    Precision volume() override;

    StaticMatrix<4, 4> integrate_NNt(const ElementBasis& basis) const;

    MapMatrix mass(Precision* buffer) override;

    Vec3 global_point_current(Precision r, Precision s) const;

    Precision integrate_scalar_field(bool scale_by_density, const ScalarField& field) override;

    Vec3 integrate_vector_field(bool scale_by_density, const VecField& field) override;

    void integrate_vector_field(Field& node_loads, bool scale_by_density, const VecField& field) override;

    Mat3 integrate_tensor_field(bool scale_by_density, const TenField& field) override;

    void compute_compliance(Field& displacement, Field& result) override;

    void compute_compliance_angle_derivative(Field& displacement,
                                             Field& result) override;

    bool compute_shear_flow(Field&       shear_flow,
                            const Field& displacement,
                            int          offset) override;

    bool compute_beam_section_forces(Field&       section_forces,
                                     const Field& displacement,
                                     int          offset) override;

    bool compute_shell_section_forces(Field&       resultants,
                                      Field&       contribution_count,
                                      const Field& displacement) override;

};

} // namespace fem::model
