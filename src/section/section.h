//
// Created by f_eggers on 21.11.2024.
//

#ifndef SECTION_H
#define SECTION_H

#include "../material/material.h"
#include "../data/region.h"

namespace fem {

struct Section {
    using Ptr = std::shared_ptr<Section>;

    material::Material::Ptr material = nullptr;
    model::ElementRegion::Ptr region = nullptr;
};


}


#endif //SECTION_H
