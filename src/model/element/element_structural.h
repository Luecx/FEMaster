/**
 * @file element_structural.h
 * @brief Declares the abstract base class for structural FEM elements.
 *
 * Structural elements extend `ElementInterface` with routines for stiffness,
 * mass, volume, load application, and stress/strain recovery.
 *
 * @see src/model/beam/beam.h
 */

#pragma once

#include "element.h"

#include "../../material/stress.h"
#include "../../material/strain.h"

namespace fem {
namespace model {

/**
 * @struct StructuralElement
 * @brief Extends `ElementInterface` with structural-specific operations.
 */
struct StructuralElement : ElementInterface {
    explicit StructuralElement(ID elem_id)
        : ElementInterface(elem_id) {}

    ~StructuralElement() override = default;

    virtual Precision volume() = 0;
    virtual MapMatrix stiffness(Precision* buffer) = 0;
    virtual MapMatrix stiffness_geom(Precision* buffer, IPData& ip_stress, int ip_start_idx) = 0;
    virtual MapMatrix mass(Precision* buffer) = 0;

    virtual Stresses stress(NodeData& displacement, std::vector<Vec3>& rst) = 0;
    virtual Strains strain(NodeData& displacement, std::vector<Vec3>& rst) = 0;
    virtual std::vector<Vec6> section_forces(NodeData& displacement) {return std::vector<Vec6>();};

    virtual void compute_stress_strain(IPData& ip_stress,
                                       IPData& ip_strain,
                                       NodeData& displacement,
                                       int ip_offset) = 0;

    virtual void compute_stress_strain_nodal(NodeData& displacement,
                                             NodeData& stress,
                                             NodeData& strain) = 0;

    virtual void apply_vload(NodeData& node_loads, Vec3 load) = 0;
    virtual void apply_tload(NodeData& node_loads, NodeData& node_temp, Precision ref_temp) = 0;
    virtual void compute_compliance(NodeData& displacement, ElementData& result) = 0;
    virtual void compute_compliance_angle_derivative(NodeData& displacement, ElementData& result) = 0;
};

} // namespace model
} // namespace fem

