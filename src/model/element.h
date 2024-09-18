#pragma once

#include "../core/core.h"
#include "../material/material.h"
#include "../math/quadrature.h"

#include <array>

namespace fem {

namespace model {

struct ElementInterface {
    const ID elem_id = 0;

    protected:
    material::Material* material = nullptr;

    public:
    ElementInterface(ID p_elem_id)
        : elem_id(p_elem_id) {}

    void set_material(material::Material* material) {
        ElementInterface::material = material;
    }

    virtual ElDofs    dofs()                                                                                     = 0;
    virtual Dim       dimensions()                                                                               = 0;
    virtual Dim       n_nodes()                                                                                  = 0;
    virtual Dim       n_integration_points()                                                                     = 0;
    virtual ID*       nodes()                                                                                    = 0;

    // iterator to iterate over nodes
    inline ID* begin() {  return nodes(); }
    inline ID* end() { return nodes() + n_nodes(); }

    virtual Precision volume(NodeData& node_coords)                                                              = 0;
    virtual MapMatrix stiffness(NodeData& position, Precision* buffer)                                           = 0;
    virtual MapMatrix mass(NodeData& position, Precision* buffer)                                                = 0;
    virtual void      compute_stress_strain_nodal(NodeData& node_coords,
                                                  NodeData& displacement,
                                                  NodeData& stress,
                                                  NodeData& strain)                                              = 0;
    virtual void      compute_stress_strain(NodeData& node_coords,
                                            NodeData& displacement,
                                            NodeData& stress,
                                            NodeData& strain,
                                            NodeData& xyz)                                                       = 0;
    virtual void      apply_dload(NodeData& node_coords, NodeData& node_loads, ID surface, StaticVector<3> load) = 0;
    virtual void      apply_vload(NodeData& node_coords, NodeData& node_loads,             StaticVector<3> load) = 0;
    virtual void      compute_compliance(NodeData& node_coords, NodeData& displacement, ElementData& result)     = 0;
};

struct ElementInterface;
using ElementPtr = std::shared_ptr<ElementInterface>;

}    // namespace model

}    // namespace fem
