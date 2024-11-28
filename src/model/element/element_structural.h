#pragma once

#include "element.h"

namespace fem::model {

struct StructuralElement : public ElementInterface {
    StructuralElement(ID p_elem_id)
        : ElementInterface(p_elem_id) {
    }
    ~StructuralElement() override = default;

    virtual Precision  volume()                                                          = 0;
    virtual MapMatrix  stiffness(Precision* buffer)                                       = 0;
    virtual MapMatrix  mass(Precision* buffer)                                            = 0;
    virtual void       compute_stress_strain_nodal(
                                                   NodeData& displacement,
                                                   NodeData& stress,
                                                   NodeData& strain)                                          = 0;
    virtual void       compute_stress_strain(
                                             NodeData& displacement,
                                             NodeData& stress,
                                             NodeData& strain,
                                             NodeData& xyz)                                                   = 0;
    virtual void       apply_vload(NodeData& node_loads, Vec3 load) = 0;
    virtual void       apply_tload(NodeData& node_loads, NodeData& node_temp, Precision ref_temp) = 0;
    virtual void       compute_compliance(NodeData& displacement, ElementData& result) = 0;
    virtual void       compute_compliance_angle_derivative(NodeData& displacement, ElementData& result) = 0;

};

}