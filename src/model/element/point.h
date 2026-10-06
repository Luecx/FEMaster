#pragma once

#include "../../section/section_point_mass.h"
#include "element_structural.h"

#include <array>
#include <memory>
#include <string>

namespace fem::model {

// One-node zero-dimensional structural element for concentrated mass, inertia
// and spring properties.
//
// PointElement stores only its element identity and one-node connectivity. All
// physical properties are supplied by an assigned PointMassSection.
//
// The element has no geometric extent, no integration points and no constitutive
// material state. Depending on its section it may contribute translational mass,
// rotary inertia, translational springs and rotational springs. Without a
// section the element is inert and activates no degrees of freedom.
//
// The spring response is linear and acts directly against ground:
//
//     f_int = K u.
//
// Consequently the complete tangent is constant, the geometric tangent is zero
// and no persistent constitutive state has to be maintained.
//
// Point elements do not provide stress, strain, PEEQ or section-resultant
// recovery. These optional capabilities use the default implementations of
// StructuralElement.
struct PointElement : StructuralElement {
    static constexpr Index N = 1;

    // Single global node connected to the point element.
    std::array<ID, N> node_ids{};

    // Construct the persistent point-element definition from element id and
    // connectivity. Section assignment, compiled offsets and ModelData binding
    // are supplied later during model compilation.
    PointElement(ID elem_id, std::array<ID, N> nodes);

    ~PointElement() override = default;

    // Create an independent copy containing only persistent element topology.
    ElementPtr copy() const override;

    // Determine active structural DOFs from mass, inertia and spring properties.
    ElDofs dofs() const override;

    // Zero-dimensional element topology.
    Dim dimensions() const override { return 0; }
    Dim n_nodes()    const override { return N; }
    Dim num_ip()     const override { return 0; }

    // Connectivity and element identification.
    const ID*   nodes()     const override;
    std::string type_name() const override;

    // A concentrated point property has no geometric volume.
    Precision volume() override;

    // Evaluate the constant spring tangent, zero geometric tangent and exact
    // internal spring force. Temperature, base state and update_state have no
    // influence because PointElement is linear and has no constitutive history.
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

    // Assemble concentrated translational mass and rotary inertia.
    MapMatrix mass(Precision* buffer) override;

    // Generic field integration uses the concentrated section mass as the
    // complete density-scaled integration measure.
    Precision integrate_scalar_field(
        bool               scale_by_density,
        const ScalarField& field
    ) override;
    Vec3 integrate_vector_field(
        bool            scale_by_density,
        const VecField& field
    ) override;
    void integrate_vector_field(
        Field&          node_loads,
        bool            scale_by_density,
        const VecField& field
    ) override;
    Mat3 integrate_tensor_field(
        bool            scale_by_density,
        const TenField& field
    ) override;
};

} // namespace fem::model