#pragma once

#include "../data/elem_data_dict.h"
#include "element.h"

namespace fem::model {

struct StructuralElement : public ElementInterface {
    StructuralElement(ID p_elem_id)
        : ElementInterface(p_elem_id) {
        _set_type(ElementTypes::StructuralType);
    }

    virtual Precision  volume(NodeData& node_coords)                                                          = 0;
    virtual MapMatrix  stiffness(NodeData& position, Precision* buffer)                                       = 0;
    virtual MapMatrix  mass(NodeData& position, Precision* buffer)                                            = 0;
    virtual void       compute_stress_strain_nodal(NodeData& node_coords,
                                                   NodeData& displacement,
                                                   NodeData& stress,
                                                   NodeData& strain)                                          = 0;
    virtual void       compute_stress_strain(NodeData& node_coords,
                                             NodeData& displacement,
                                             NodeData& stress,
                                             NodeData& strain,
                                             NodeData& xyz)                                                   = 0;
    virtual void       apply_vload(NodeData& node_coords, NodeData& node_loads, Vec3 load) = 0;
    virtual void       apply_tload(NodeData& node_coords, NodeData& node_loads, NodeData& node_temp, Precision ref_temp) = 0;
    virtual void       compute_compliance(NodeData& node_coords, NodeData& displacement, ElementData& result) = 0;
};

}