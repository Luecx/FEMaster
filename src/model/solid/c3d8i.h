/**
 * @file c3d8i.h
 * @brief Eight-node fully integrated hexahedron with incompatible modes.
 *
 * C3D8I augments the compatible trilinear deformation gradient with thirteen
 * element-local enhanced modes. Nine vectorial modes improve bending/shear
 * response and four scalar volumetric modes reduce near-incompressible locking.
 * The internal variables are solved locally and statically condensed.
 *
 * Finite deformation uses an incremental enhanced deformation gradient updated
 * multiplicatively from the last accepted configuration. The enhanced history is
 * committed only with an accepted nonlinear load increment.
 */

#pragma once

#include "c3d8.h"

#include <array>
#include <vector>

namespace fem::model {

class C3D8I final : public C3D8 {
public:
    static constexpr Index N          = 8;
    static constexpr Dim   D          = 3;
    static constexpr Index ndof       = N * D;
    static constexpr Index n_internal = 13;
    static constexpr Index n_total    = ndof + n_internal;

    using Vector13   = StaticVector<n_internal>;
    using Vector24   = StaticVector<ndof>;
    using Vector37   = StaticVector<n_total>;
    using Matrix13   = StaticMatrix<n_internal, n_internal>;
    using Matrix24   = StaticMatrix<ndof, ndof>;
    using Matrix37   = StaticMatrix<n_total, n_total>;
    using Matrix6x13 = StaticMatrix<6, n_internal>;
    using Matrix6x37 = StaticMatrix<6, n_total>;

    C3D8I(ID elem_id, const std::array<ID, N>& node_ids);
    ~C3D8I() override = default;

    ElementPtr copy() const override { return std::make_shared<C3D8I>(elem_id, node_ids); }
    std::string type_name() const override { return "C3D8I"; }

    MapMatrix stiffness(Precision* buffer) override;
    MapMatrix stiffness_geom(Precision* buffer, const Field& displacement) override;
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
    struct TrialState {
        StaticMatrix<N, D> coordinates = StaticMatrix<N, D>::Zero();
        std::array<Mat3, N> deformation_gradients{};
        bool valid = false;
    };

    struct LinearData {
        Matrix24 Kuu = Matrix24::Zero();
        StaticMatrix<ndof, n_internal> Kua =
            StaticMatrix<ndof, n_internal>::Zero();
        StaticMatrix<n_internal, ndof> Kau =
            StaticMatrix<n_internal, ndof>::Zero();
        Matrix13 Kaa = Matrix13::Zero();
        std::array<StaticMatrix<6, ndof>, N> B{};
        std::array<Matrix6x13, N> M{};
        Precision characteristic_length = Precision(1);
    };

    struct NonlinearData {
        Vector37 residual = Vector37::Zero();
        Matrix37 tangent  = Matrix37::Zero();
        std::array<Mat3, N> deformation_gradients{};
    };

    LinearData linear_data();
    Vector13 linear_internal_parameters(const LinearData& data, const Field& displacement);

    Matrix6x13 incompatible_strain_modes(
        const StaticMatrix<N, D>& coordinates,
        Precision                 r,
        Precision                 s,
        Precision                 t,
        Precision                 characteristic_length
    );

    NonlinearData evaluate_nonlinear(
        const StaticMatrix<N, D>& current_coordinates,
        const Vector13&           internal,
        bool                      write_material_state
    );

    NonlinearData solve_internal_modes(
        const StaticMatrix<N, D>& current_coordinates,
        Vector13&                 internal
    );

    void write_material_state(const std::array<Mat3, N>& deformation_gradients);
    void scatter_force(NodeData& nodal_forces, const Vector24& local_force);

    static Vec6 strain_variation(const Mat3& F, const Mat3& dF);
    static Precision stress_second_variation(
        const Mat3& S,
        const Mat3& F,
        const Mat3& dFp,
        const Mat3& dFq,
        const Mat3& d2F
    );

    Precision reference_characteristic_length();

    bool initialized_ = false;
    StaticMatrix<N, D> committed_coordinates_ = StaticMatrix<N, D>::Zero();
    std::array<Mat3, N> committed_deformation_gradients_{};
    std::vector<TrialState> trial_stack_;
};

} // namespace fem::model
