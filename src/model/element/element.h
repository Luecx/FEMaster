#pragma once

#include "../../core/core.h"
#include "../../material/material.h"
#include "../geometry/surface/surface.h"
#include "../../math/quadrature.h"
#include "element_types.h"

#include "../data/node_data_dict.h"
#include "../data/elem_data_dict.h"

#include <array>

namespace fem {

namespace model {

struct ElementInterface : public std::enable_shared_from_this<ElementInterface> {
    const ID elem_id = 0;

    protected:
    // storing the type of this element
    ElementTypeFlags _flags;

    material::Material::Ptr _material = nullptr;
    NodeDataDict::Ptr       _node_data_dict = nullptr;
    ElemDataDict::Ptr       _elem_data_dict = nullptr;

    public:
    ElementInterface(ID p_elem_id)
        : elem_id(p_elem_id) {}

    virtual ~ElementInterface() = default;

    void set_material(material::Material::Ptr material) {
        _material = material;
    }
    void set_node_data_dict(NodeDataDict::Ptr node_data) {
        _node_data_dict = node_data;
    }
    void set_elem_data_dict(ElemDataDict::Ptr elem_data) {
        _elem_data_dict = elem_data;
    }

    virtual ElDofs dofs()                 = 0;
    virtual Dim    dimensions()           = 0;
    virtual Dim    n_nodes()              = 0;
    virtual Dim    n_integration_points() = 0;
    virtual ID*    nodes()                = 0;

    virtual SurfacePtr surface(ID surface_id) = 0;

    // iterator to iterate over nodes
    inline ID* begin() {
        return nodes();
    }
    inline ID* end() {
        return nodes() + n_nodes();
    }

    // cast to a specific type
    // Templated function for casting to a specific type
    template <typename T>
    std::shared_ptr<T> as() {
        return std::dynamic_pointer_cast<T>(shared_from_this());
    }

    void _set_type(ElementType flags) {
        _flags |= flags;
    }
    bool is_type(ElementType type) const {
        return (_flags & type) != 0;
    }

};

struct ElementInterface;
using ElementPtr = std::unique_ptr<ElementInterface>;

}    // namespace model

}    // namespace fem
