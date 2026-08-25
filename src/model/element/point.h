/**
 * @file point.h
 * @brief Declares a one-node zero-dimensional structural point element.
 *
 * `PointElement` owns only element identity and one-node connectivity. Physical
 * point properties are supplied by `PointMassSection`, exactly like other
 * structural elements obtain their properties from their assigned section.
 *
 * Without a point-mass section the element is inert and activates no DOFs. With
 * a section it contributes concentrated translational mass, rotary inertia and
 * diagonal translational/rotational stiffness against ground. Density-scaled
 * field integration interprets the section mass as the complete concentrated
 * mass measure, which makes element-based gravity and inertia loads work through
 * the ordinary structural-element path.
 *
 * @see StructuralElement
 * @see PointMassSection
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#pragma once

#include "element_structural.h"
#include "../../section/section_point_mass.h"

#include <array>
#include <memory>
#include <string>

namespace fem::model {

/**
 * @struct PointElement
 * @brief One-node structural element using a concentrated point-mass section.
 *
 * The element has no geometric integration measure and no constitutive state.
 * Its section directly defines diagonal mass/inertia and ground-stiffness
 * coefficients. Missing sections are permitted and leave the topology entirely
 * inactive.
 */
struct PointElement : StructuralElement {
    static constexpr Index N = 1;

    std::array<ID, N> node_ids {};

    PointElement(ID elem_id, std::array<ID, N> nodes)
        : StructuralElement(elem_id), node_ids(nodes) {}

    ~PointElement() override = default;

    ElementPtr copy() const override {
        return std::make_shared<PointElement>(elem_id, node_ids);
    }

    ElDofs dofs() const override {
        if (!_section) return ElDofs{false, false, false, false, false, false};

        const auto* section = _section->as<PointMassSection>();
        logging::error(section != nullptr,
            "PointElement: section is not a PointMassSection for element ", elem_id);

        return ElDofs{
            section->mass_ != Precision(0) || section->spring_constants_(0) != Precision(0),
            section->mass_ != Precision(0) || section->spring_constants_(1) != Precision(0),
            section->mass_ != Precision(0) || section->spring_constants_(2) != Precision(0),
            section->rotary_inertia_(0) != Precision(0) || section->rotary_spring_constants_(0) != Precision(0),
            section->rotary_inertia_(1) != Precision(0) || section->rotary_spring_constants_(1) != Precision(0),
            section->rotary_inertia_(2) != Precision(0) || section->rotary_spring_constants_(2) != Precision(0)
        };
    }

    Dim dimensions() const override { return 0; }
    Dim n_nodes() const override { return N; }
    Dim num_ip() const override { return 0; }
    const ID* nodes() const override { return node_ids.data(); }
    std::string type_name() const override { return "POINT"; }

    Precision volume() override { return Precision(0); }

    MapMatrix stiffness(Precision* buffer) override {
        MapMatrix result(buffer, 6, 6);
        result.setZero();
        if (!_section) return result;

        const auto* section = _section->as<PointMassSection>();
        logging::error(section != nullptr,
            "PointElement: section is not a PointMassSection for element ", elem_id);

        result(0, 0) = section->spring_constants_(0);
        result(1, 1) = section->spring_constants_(1);
        result(2, 2) = section->spring_constants_(2);
        result(3, 3) = section->rotary_spring_constants_(0);
        result(4, 4) = section->rotary_spring_constants_(1);
        result(5, 5) = section->rotary_spring_constants_(2);
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
        internal_force_nonlinear(ip_stress_state, nodal_forces, displacement);
        return stiffness(buffer);
    }

    MapMatrix mass(Precision* buffer) override {
        MapMatrix result(buffer, 6, 6);
        result.setZero();
        if (!_section) return result;

        const auto* section = _section->as<PointMassSection>();
        logging::error(section != nullptr,
            "PointElement: section is not a PointMassSection for element ", elem_id);

        result(0, 0) = section->mass_;
        result(1, 1) = section->mass_;
        result(2, 2) = section->mass_;
        result(3, 3) = section->rotary_inertia_(0);
        result(4, 4) = section->rotary_inertia_(1);
        result(5, 5) = section->rotary_inertia_(2);
        return result;
    }

    void internal_force_nonlinear(Field& ip_stress_state,
                                  NodeData& nodal_forces,
                                  const Field& displacement) override {
        (void) ip_stress_state;
        if (!_section) return;

        const auto* section = _section->as<PointMassSection>();
        logging::error(section != nullptr,
            "PointElement: section is not a PointMassSection for element ", elem_id);

        const Index node = static_cast<Index>(node_ids[0]);
        for (Index dof = 0; dof < 3; ++dof) {
            nodal_forces(node, dof)     += section->spring_constants_(dof) * displacement(node, dof);
            nodal_forces(node, dof + 3) += section->rotary_spring_constants_(dof) * displacement(node, dof + 3);
        }
    }

    Precision integrate_scalar_field(bool scale_by_density, const ScalarField& field) override {
        if (!scale_by_density || !_section) return Precision(0);
        const auto* section = _section->as<PointMassSection>();
        logging::error(section != nullptr,
            "PointElement: section is not a PointMassSection for element ", elem_id);
        return section->mass_ * field(node_position(0));
    }

    Vec3 integrate_vector_field(bool scale_by_density, const VecField& field) override {
        if (!scale_by_density || !_section) return Vec3::Zero();
        const auto* section = _section->as<PointMassSection>();
        logging::error(section != nullptr,
            "PointElement: section is not a PointMassSection for element ", elem_id);
        return section->mass_ * field(node_position(0));
    }

    void integrate_vector_field(Field& node_loads, bool scale_by_density, const VecField& field) override {
        if (!scale_by_density || !_section) return;
        const auto* section = _section->as<PointMassSection>();
        logging::error(section != nullptr,
            "PointElement: section is not a PointMassSection for element ", elem_id);

        const Index node = static_cast<Index>(node_ids[0]);
        const Vec3 value = section->mass_ * field(node_position(0));
        node_loads(node, 0) += value(0);
        node_loads(node, 1) += value(1);
        node_loads(node, 2) += value(2);
    }

    Mat3 integrate_tensor_field(bool scale_by_density, const TenField& field) override {
        if (!scale_by_density || !_section) return Mat3::Zero();
        const auto* section = _section->as<PointMassSection>();
        logging::error(section != nullptr,
            "PointElement: section is not a PointMassSection for element ", elem_id);
        return section->mass_ * field(node_position(0));
    }

    void apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) override {
        (void) node_loads;
        (void) node_temp;
        (void) ref_temp;
    }

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
