//
// Created by f_eggers on 21.11.2024.
//

#ifndef SECTION_H
#define SECTION_H

#include "../material/material.h"
#include "../data/region.h"

#include "profile.h"

namespace fem {

struct Section {
    using Ptr = std::shared_ptr<Section>;

    virtual ~Section() = default; // Make Section polymorphic

    material::Material::Ptr material = nullptr;
    model::ElementRegion::Ptr region = nullptr;

    template<typename T>
    T* as() {
        return dynamic_cast<T*>(this);
    }
};

struct BeamSection : Section{
    using Ptr = std::shared_ptr<BeamSection>;
    Vec3 n1;
    Profile::Ptr profile = nullptr;
};


}


#endif //SECTION_H
