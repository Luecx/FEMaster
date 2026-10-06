/**
 * @file b31.h
 * @brief Defines the two-node geometrically exact three-dimensional beam.
 *
 * B31 is a shear-flexible finite-rotation beam with six generalized section
 * strains. Its mechanical state is evaluated from the undeformed reference
 * geometry and the supplied total nodal displacement/rotation vectors.
 */

#pragma once

#include "beam.h"
#include "../geometry/line/line2a.h"

#include <array>

namespace fem::model {

/**
 * @brief Two-node geometrically exact shear-flexible beam.
 *
 * The element uses two endpoint (GLL order-one) integration points. Nodal
 * rotations follow FEMaster's global total axis-angle convention. The section
 * reference line is mapped to the shear point with an exact rotating rigid
 * offset, while the section constitutive matrix retains the centroid/shear-point
 * coupling defined by Profile.
 */
struct B31 : BeamElement<2> {
    static constexpr Index N             = 2;
    static constexpr Index dofs_per_node = 6;
    static constexpr Index num_dofs      = N * dofs_per_node;
    static constexpr Index n_strain      = 6;

    using Vec12  = StaticVector<num_dofs>;
    using Mat12  = StaticMatrix<num_dofs, num_dofs>;
    using Mat6x12 = StaticMatrix<n_strain, num_dofs>;

    B31(ID elem_id, std::array<ID, N> node_ids_in)
        : BeamElement(elem_id, node_ids_in) {}

    ElementPtr copy() const override {
        return std::make_shared<B31>(elem_id, node_ids);
    }

    std::string type_name() const override { return "B31"; }

    StaticMatrix<12, 12> stiffness_impl() override;
    StaticMatrix<12, 12> stiffness_geom_impl(
        const Field& target_displacement,
        const Field* target_temperature,
        const Field* base_temperature
    ) override;
    StaticMatrix<12, 12> mass_impl() override;

    MapMatrix evaluate(
        Precision*   tangent,
        Precision*   geometric_tangent,
        NodeData*    internal_force,
        const Field* target_displacement,
        const Field* target_temperature,
        const Field* base_displacement,
        const Field* base_temperature,
        bool         update_state
    ) override;

    Dim num_ip() const override { return 2; }

    RowMatrix stress_strain_nodal_rst() override;
    RowMatrix stress_strain_ip_rst() override;

    void compute_stress_strain(
        Field*           strain,
        Field*           stress,
        const Field&     target_displacement,
        const Field*     target_temperature,
        const RowMatrix& rst,
        const Field*     base_displacement,
        const Field*     base_temperature
    ) override;

    bool compute_beam_section_forces(
        Field&       section_forces,
        const Field& target_displacement,
        const Field* target_temperature,
        const Field* base_displacement = nullptr,
        const Field* base_temperature = nullptr
    ) override;

    LinePtr line(ID line_id) override {
        (void) line_id;
        return std::make_shared<Line2A>(this->node_ids);
    }

private:
    struct PointKinematics {
        Vec6      strain = Vec6::Zero();
        Mat6x12   B      = Mat6x12::Zero();
    };

    struct Kinematics {
        std::array<PointKinematics, N> point;
    };

    Vec3 node_position_reference(Index local_node) const;
    Precision reference_length() const;
    Mat3 reference_frame();
    Vec12 element_displacement(const Field* displacement) const;

    Kinematics kinematics(const Vec12& q, bool with_B);
    StaticMatrix<6, 6> constitutive_matrix();
    Vec6 thermal_strain(Index point, const Field* temperature);
    std::array<Vec6, N> section_resultants(
        const Kinematics& state,
        const Field* temperature
    );

    Vec12 exact_internal_force(
        const Vec12& q,
        const Field* temperature
    );

    Mat12 material_tangent(const Kinematics& state);
    Mat12 geometric_tangent_from_resultants(
        const Vec12& q,
        const std::array<Vec6, N>& resultants
    );

    std::array<Vec6, N> continued_resultants(
        const Kinematics& base_state,
        const std::array<Vec6, N>& base_resultants,
        const std::array<Vec6, N>& target_temperature_resultants,
        const Vec12& delta
    );

    Vec6 resultant_about_reference(const Vec6& shear_point_resultant) const;
};

} // namespace fem::model
