/**
 * @file point.h
 * @brief Declares a one-node zero-dimensional structural point element.
 *
 * `PointElement` preserves ordinary element identity, ELSET membership and
 * Part/Instance expansion for zero-dimensional element topologies. The element
 * itself carries no physical property: stiffness, mass, constitutive response
 * and density-scaled integration all evaluate to zero.
 *
 * Abaqus `MASS` properties use the point element only to resolve target nodes and
 * create separate `feature::PointMass` objects. Keeping the property outside this
 * class avoids turning a topology carrier into a second point-mass implementation
 * while still retaining real MASS element identifiers and sets.
 *
 * @see StructuralElement
 * @see feature::PointMass
 *
 * @author Finn Eggers
 * @date 24.08.2026
 */

#pragma once

#include "element_structural.h"

#include <array>
#include <memory>
#include <string>

namespace fem::model {

/**
 * @struct PointElement
 * @brief One-node structural topology carrier without direct physics.
 *
 * The class deliberately implements the complete `StructuralElement` interface
 * with zero or no-op behavior. It therefore participates in the regular element
 * lifecycle and compiled topology without owning a section, material, mass or
 * integration scheme. Physical point properties remain separate model features.
 */
struct PointElement : StructuralElement {
    static constexpr Index N = 1;

    std::array<ID, N> node_ids {};

    PointElement(ID elem_id, std::array<ID, N> nodes)
        : StructuralElement(elem_id), node_ids(nodes) {}

    ~PointElement() override = default;

    // Preserve only element identity and connectivity through Instance expansion.
    ElementPtr copy() const override {
        return std::make_shared<PointElement>(elem_id, node_ids);
    }

    // The point topology itself activates no generalized DOFs.
    ElDofs    dofs() const override { return ElDofs{false, false, false, false, false, false}; }
    Dim       dimensions() const override { return 0; }
    Dim       n_nodes() const override { return N; }
    Dim       num_ip() const override { return 0; }
    const ID* nodes() const override { return node_ids.data(); }
    std::string type_name() const override { return "POINT"; }

    // A zero-dimensional topology has no geometric measure or element operators.
    Precision volume() override { return Precision(0); }

    MapMatrix stiffness(Precision* buffer) override {
        MapMatrix result(buffer, 6, 6);
        result.setZero();
        return result;
    }

    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override {
        (void) ip_stress;
        (void) ip_start_idx;
        MapMatrix result(buffer, 6, 6);
        result.setZero();
        return result;
    }

    MapMatrix stiffness_tangent(Precision* buffer,
                                Field& ip_stress_state,
                                NodeData& nodal_forces,
                                const Field& displacement) override {
        (void) ip_stress_state;
        (void) nodal_forces;
        (void) displacement;
        return stiffness(buffer);
    }

    MapMatrix mass(Precision* buffer) override {
        MapMatrix result(buffer, 6, 6);
        result.setZero();
        return result;
    }

    void internal_force_nonlinear(Field& ip_stress_state,
                                  NodeData& nodal_forces,
                                  const Field& displacement) override {
        (void) ip_stress_state;
        (void) nodal_forces;
        (void) displacement;
    }

    // Point topology has no integration measure. Concentrated mass integration
    // is handled by the separately created PointMass feature.
    Precision integrate_scalar_field(bool scale_by_density, const ScalarField& field) override {
        (void) scale_by_density;
        (void) field;
        return Precision(0);
    }

    Vec3 integrate_vector_field(bool scale_by_density, const VecField& field) override {
        (void) scale_by_density;
        (void) field;
        return Vec3::Zero();
    }

    void integrate_vector_field(Field& node_loads, bool scale_by_density, const VecField& field) override {
        (void) node_loads;
        (void) scale_by_density;
        (void) field;
    }

    Mat3 integrate_tensor_field(bool scale_by_density, const TenField& field) override {
        (void) scale_by_density;
        (void) field;
        return Mat3::Zero();
    }

    void apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) override {
        (void) node_loads;
        (void) node_temp;
        (void) ref_temp;
    }

    // Point elements have no strain, stress or internal-force state.
    void compute_stress_strain(Field* strain,
                               Field* stress,
                               const Field& displacement,
                               const RowMatrix& rst,
                               int offset,
                               bool use_green_lagrange_nl) override {
        (void) strain;
        (void) stress;
        (void) displacement;
        (void) rst;
        (void) offset;
        (void) use_green_lagrange_nl;
    }

    void compute_stress_state(Field& stress_state,
                              const Field& displacement,
                              int offset,
                              bool use_green_lagrange_nl) override {
        (void) stress_state;
        (void) displacement;
        (void) offset;
        (void) use_green_lagrange_nl;
    }

    void compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) override {
        (void) node_forces;
        (void) ip_stress;
    }
};

} // namespace fem::model
