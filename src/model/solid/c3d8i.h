/**
 * @file c3d8i.h
 * @brief Eight-node incompatible-mode hexahedral solid.
 *
 * C3D8I augments the full-integration trilinear brick by thirteen element-local
 * incompatible deformation modes. Nine Wilson modes relax bending and four
 * scalar volumetric modes reduce near-incompressible locking. The internal
 * variables are solved locally and statically condensed from the global tangent.
 *
 * Finite deformation uses an incremental enhanced deformation gradient. The
 * enhanced gradient is committed only with an accepted nonlinear load increment,
 * so every Newton and line-search evaluation starts from the same physical
 * beginning-of-increment configuration.
 */

#pragma once

#include "c3d8.h"

#include <array>

namespace fem::model {

class C3D8I final : public C3D8 {
public:
    static constexpr Index N                 = 8;
    static constexpr Dim   D                 = 3;
    static constexpr Index ndof              = 24;
    static constexpr Index principal_modes   = 9;
    static constexpr Index volumetric_modes  = 4;
    static constexpr Index internal_dofs     = 13;
    static constexpr Index total_dofs        = ndof + internal_dofs;

    using Matrix24 = StaticMatrix<ndof, ndof>;
    using Vector24 = StaticVector<ndof>;
    using Matrix13 = StaticMatrix<internal_dofs, internal_dofs>;
    using Vector13 = StaticVector<internal_dofs>;
    using Matrix37 = StaticMatrix<total_dofs, total_dofs>;
    using Vector37 = StaticVector<total_dofs>;
    using BMatrix  = StaticMatrix<6, total_dofs>;

    C3D8I(ID elem_id, const std::array<ID, N>& node_ids);
    ~C3D8I() override = default;

    ElementPtr copy() const override { return std::make_shared<C3D8I>(elem_id, node_ids); }
    std::string type_name() const override { return "C3D8I"; }

    MapMatrix stiffness(Precision* buffer) override;
    MapMatrix stiffness_tangent(
        Precision*   buffer,
        NodeData&    nodal_forces,
        const Field& displacement
    ) override;

    void compute_stress_strain(
        Field*           strain,
        Field*           stress,
        const Field&     displacement,
        const RowMatrix& rst,
        int              offset,
        bool             use_green_lagrange_nl
    ) override;

    void step_begin() override;
    void step_end() override;
    void nonlinear_begin_increment() override;
    void nonlinear_commit_increment() override;
    void nonlinear_rollback_increment() override;

private:
    struct NonlinearEvaluation {
        Vector37 residual = Vector37::Zero();
        Matrix37 tangent  = Matrix37::Zero();
        std::array<Mat3, N> enhanced_F {};
    };

    struct InternalSolution {
        Vector13 parameters = Vector13::Zero();
        NonlinearEvaluation evaluation;
    };

    Precision characteristic_length_ = Precision(1);
    StaticMatrix<N, D> committed_coords_ = StaticMatrix<N, D>::Zero();
    StaticMatrix<N, D> trial_coords_     = StaticMatrix<N, D>::Zero();
    std::array<Mat3, N> committed_F_ {};
    std::array<Mat3, N> trial_F_ {};
    Vector13 trial_internal_ = Vector13::Zero();
    bool nonlinear_state_initialized_ = false;

    static Vec6 strain_variation(const Mat3& F, const Mat3& dF);
    static Vec6 linearized_strain(const Mat3& gradient);
    static Precision stress_contraction(const Mat3& stress, const Mat3& variation);

    StaticMatrix<6, internal_dofs> linear_internal_B(
        Precision r,
        Precision s,
        Precision t,
        const StaticMatrix<N, D>& reference_coords,
        Precision& det0
    );

    Matrix37 linear_full_stiffness();
    Vector13 linear_internal_parameters(const Field& displacement);
    NonlinearEvaluation evaluate_nonlinear(
        const Field&    displacement,
        const Vector13& internal,
        bool            with_tangent,
        bool            write_material_state
    );
    InternalSolution solve_internal(
        const Field& displacement,
        bool         with_tangent,
        bool         write_material_state,
        bool         update_trial_cache
    );

    void scatter_force(NodeData& nodal_forces, const Vector24& local_force);
};

} // namespace fem::model
