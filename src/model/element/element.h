#pragma once

#include "../../core/core.h"
#include "../../material/material.h"
#include "../geometry/surface/surface.h"
#include "../../math/quadrature.h"
#include "element_types.h"

#include <array>

namespace fem {

namespace model {

struct ElementInterface : public std::enable_shared_from_this<ElementInterface> {
    const ID elem_id = 0;

    protected:
    // storing the type of this element
    ElementTypeFlags _flags;

    material::Material* _material = nullptr;

    public:
    ElementInterface(ID p_elem_id)
        : elem_id(p_elem_id) {}

    virtual ~ElementInterface() = default;


    void set_material(material::Material* material) {
        ElementInterface::_material = material;
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

    // Function to check if the element is of a specific type
    bool is_type(ElementType type) const {
        return (_flags & type) != 0;
    }

};

struct ElementInterface;
using ElementPtr = std::shared_ptr<ElementInterface>;

}    // namespace model

}    // namespace fem
