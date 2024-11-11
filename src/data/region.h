//
// Created by Finn Eggers on 02.11.24.
//

#ifndef FEMASTER_REGION_H
#define FEMASTER_REGION_H

#include <memory>

#include "region_type.h"
#include "collection.h"
#include "../core/core.h"

namespace fem::model {

template<RegionTypes RT>
struct Region : public Collection<ID> {
    using Ptr = std::shared_ptr<Region<RT>>;

    Region(std::string name) : Collection(name, true, false) {}
};

using NodeRegion    = Region<RegionTypes::NODE>;
using ElementRegion = Region<RegionTypes::ELEMENT>;
using SurfaceRegion = Region<RegionTypes::SURFACE>;

}

#endif    // FEMASTER_REGION_H
