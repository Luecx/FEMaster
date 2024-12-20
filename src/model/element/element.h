#pragma once

#include "../model_data.h"
#include "../geometry/surface/surface.h"
#include "../../section/section.h"

namespace fem {

namespace model {

struct ElementInterface {
    const ID elem_id = 0;

    Section::Ptr _section = nullptr;
    ModelDataPtr _model_data;

    ElementInterface(ID p_elem_id)
        : elem_id(p_elem_id) {}

    virtual ~ElementInterface() {};


    void set_section(Section::Ptr section) {
        _section = section;
    }
    material::MaterialPtr material() {
        return _section->material;
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
    // Templated function for casting to a specific type
    template <typename T>
    T* as() {
        return dynamic_cast<T*>(this);
    }
};

}    // namespace model

}    // namespace fem
