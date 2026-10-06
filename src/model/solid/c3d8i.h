/**
 * @file c3d8i.h
 * @brief Declares the incompatible-mode eight-node hexahedral solid.
 *
 * C3D8I extends the fully integrated C3D8 continuum element by thirteen
 * element-local enhanced deformation modes. The compatible trilinear geometry,
 * topology, faces, mass integration and common solid utilities remain inherited
 * from C3D8.
 *
 * The enhanced parameters never become global degrees of freedom. Every
 * expansion point uses right multiplicative reference enhancement, a local
 * stationarity solve and Schur-complement condensation of the coupled tangent.
 * Linear mechanics evaluates the same formulation at zero nodal displacement.
 *
 * Temperature, stress recovery and perturbation geometric stiffness use the
 * identical enhanced base state so all element operators represent the same
 * thermo-mechanical kinematic approximation.
 *
 * @see C3D8
 * @see SolidElement
 *
 * @author Finn Eggers
 * @date 02.10.2026
 */

#pragma once

#include "c3d8.h"

#include <array>

namespace fem::model {

/**
 * @brief Fully integrated eight-node solid with thirteen incompatible modes.
 *
 * C3D8I retains the external C3D8 topology and its 24 translational nodal
 * degrees of freedom while enriching the element kinematics by thirteen local
 * deformation-gradient modes. Nine vector-valued principal modes improve the
 * low-order bending response and four scalar volumetric modes reduce excessive
 * volumetric constraint.
 *
 * The enhanced parameters are element-local unknowns. The tangent equations at
 * each expansion point have the block form
 *
 * [ Kuu  Kua ] [ u     ] = [ fu ]
 * [ Kau  Kaa ] [ alpha ]   [ fa ],
 *
 * and enhanced increments are eliminated locally. Every evaluation uses
 * the objective kinematics
 *
 * F_bar = F_c (I + sum_m alpha_m H_m)
 *
 * and the same block structure for the consistent Total-Lagrangian tangent
 * after a local Newton solve of the enhanced residual r_alpha = 0. The enhanced
 * tensors H_m remain in the reference configuration, so a superposed spatial
 * rotation maps F_bar directly to Q F_bar.
 *
 * Material history follows the common structural-element contract. Auxiliary
 * stiffness, perturbation and recovery paths read only committed material state.
 * Only evaluate() writes the converged constitutive trial state.
 *
 * No enhanced parameter is stored persistently in the element object; every
 * evaluation reconstructs the local stationary state from the supplied nodal
 * configuration.
 */
class C3D8I final : public C3D8 {
public:
    // Fixed topology and local enhanced-mode dimensions
    static constexpr Index N       = 8;
    static constexpr Dim   D       = 3;
    static constexpr Index ndof    = N * D;
    static constexpr Index n_modes = 13;

    // Fixed-size algebra used by the local EAS equations and condensation
    using Vector24      = StaticVector<ndof>;
    using Vector13      = StaticVector<n_modes>;
    using Matrix24      = StaticMatrix<ndof, ndof>;
    using Matrix13      = StaticMatrix<n_modes, n_modes>;
    using Matrix24x13   = StaticMatrix<ndof, n_modes>;
    using Matrix13x24   = StaticMatrix<n_modes, ndof>;
    using Matrix6x13    = StaticMatrix<6, n_modes>;
    using EnhancedModes = std::array<Mat3, n_modes>;

    /**
     * @brief Local nodal/enhanced residual and tangent block system.
     *
     * The structure represents the element equations before static
     * condensation. kuu, kua, kau and kaa are the tangent blocks with respect
     * to nodal displacement u and enhanced parameters alpha; ru and ra are
     * their corresponding residual vectors.
     *
     * The finite formulation fills the requested residual and tangent blocks,
     * including material and geometric derivatives. Its reference evaluation
     * also supplies linear stiffness through condensation of these same blocks.
     *
     * The structure owns no persistent element state. It exists only during one
     * element evaluation and is discarded after condensation.
     */
    struct EnhancedSystem {
        // Coupled nodal/enhanced tangent blocks before static condensation
        Matrix24    kuu = Matrix24::Zero();
        Matrix24x13 kua = Matrix24x13::Zero();
        Matrix13x24 kau = Matrix13x24::Zero();
        Matrix13    kaa = Matrix13::Zero();

        // Nodal and enhanced residuals evaluated at the current local candidate
        Vector24    ru  = Vector24::Zero();
        Vector13    ra  = Vector13::Zero();
    };

    /**
     * @brief Fixed geometry of one constitutive point during a local EAS solve.
     *
     * The reference derivatives, physical reference-volume weight and enhanced
     * basis are independent of the local parameters. The compatible deformation
     * gradient is also fixed while Newton solves the enhanced stationarity
     * equations at a supplied nodal configuration. These quantities are built
     * once per element evaluation and reused by local Newton, final assembly
     * and nonlinear recovery without storing trial history in the element.
     */
    struct NonlinearPoint {
        // Natural coordinates and complete reference-volume quadrature weight
        Vec3      natural = Vec3::Zero();
        Precision measure = Precision(0);

        // Reference interpolation derivatives and compatible finite deformation
        StaticMatrix<N, D> derivatives = StaticMatrix<N, D>::Zero();
        Mat3               compatible  = Mat3::Identity();

        // Thirteen reference-configuration enhanced gradient tensors
        EnhancedModes modes;
    };

    // Full C3D8 quadrature geometry owned by one evaluation, not by the element
    using NonlinearPoints = std::array<NonlinearPoint, 8>;

    // Construction and element identity
    C3D8I(ID elem_id, const std::array<ID, N>& node_ids);
    ~C3D8I() override = default;

    ElementPtr copy() const override { return std::make_shared<C3D8I>(elem_id, node_ids); }
    std::string type_name() const override;

    // Common mechanical response with local elimination of all thirteen
    // enhanced parameters before exposing nodal operators to the global solver.
    MapMatrix evaluate(
        Precision*   tangent,
        Precision*   geometric_tangent,
        NodeData*    internal_force,
        const Field* displacement,
        const Field* linearization,
        bool         update_state
    ) override;

    // Compliance orientation sensitivity uses the stationary enhanced strain
    // rather than the compatible C3D8 strain inherited by the common solid.
    void compute_compliance_angle_derivative(Field& displacement, Field& result) override;

    // Stress and strain recovery reconstructs the stationary enhanced state
    // before evaluating the constitutive response at the requested locations.
    void compute_stress_strain(
        Field*           strain,
        Field*           stress,
        const Field&     displacement,
        const RowMatrix& rst,
        int              offset,
        const Field*     linearization) override;

private:
    // Condensed geometric stiffness generated by the linearized stress increment
    // from the base state u0 to the requested state u.
    MapMatrix stiffness_geom(
        Precision*              buffer,
        const NonlinearPoints&  points,
        const Vector13&          alpha,
        const EnhancedSystem&    system,
        const Vector24&          displacement_increment,
        const StaticVector<N>*   thermal_strain);

    // Enhanced deformation-gradient basis and its linearized or finite-strain
    // work-conjugate strain matrices.
    EnhancedModes enhanced_gradient_modes(
        const StaticMatrix<N, D>& reference_coords,
        Precision                 r,
        Precision                 s,
        Precision                 t
    );
    Matrix6x13 enhanced_green_lagrange_matrix(const Mat3& deformation_gradient, const EnhancedModes& modes);

    // Finite block-system assembly at a stationary reference or supplied state.
    // Requested residuals and tangents use the same enhanced kinematics.
    EnhancedSystem assemble_reference_system();
    EnhancedSystem assemble_nonlinear_system(
        const NonlinearPoints& points,
        const Vector13&        alpha,
        bool                   write_material_state,
        bool                   assemble_global_blocks,
        bool                   assemble_tangent,
        bool                   include_geometric = true,
        const StaticVector<N>* thermal_strain = nullptr
    );
    NonlinearPoints nonlinear_points(
        const StaticMatrix<N, D>& reference_coords,
        const StaticMatrix<N, D>& current_coords
    );

    // Solve finite local stationarity or its reference tangent increment without
    // introducing any global degrees of freedom.
    Vector13 solve_linear_modes(const Vector24& displacement);
    Vector13 solve_nonlinear_modes(
        const NonlinearPoints& points,
        const StaticVector<N>* thermal_strain = nullptr);

    // Current nodal free thermal strain used by the finite constitutive state.
    StaticVector<N> nodal_thermal_strain();

    // Map between the common global nodal fields and the element-local
    // translational ordering [u1x,u1y,u1z,...,u8x,u8y,u8z].
    Vector24 local_displacement(const Field& displacement);
    void assemble_local_force(Field& node_forces, const Vector24& local_force);
};

} // namespace fem::model
