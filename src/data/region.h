//
// Created by Finn Eggers on 02.11.24.
//

#ifndef FEMASTER_REGION_H
#define FEMASTER_REGION_H

#include <memory>

#include "../core/logging.h"
#include "../core/types_num.h"

#include "region_type.h"
#include "collection.h"

namespace fem::model {

template<RegionTypes RT>
struct Region : public Collection<ID> {
    using Ptr = std::shared_ptr<Region<RT>>;

    Region(std::string name)
        : Collection(name, true, false) {}

    void info();

    // string output operator
    std::ostream& operator<<(std::ostream& os) const {
        os << "Region: " << this->name;
        os << "   Type: " << RT;
        os << "   Size: " << this->size();
        os << "   IDs : ";
        // for(size_t i = 0; i < std::max((size_t)4, this->size()); i++) {
        //     os << this->at(i) << " ";
        // }
        return os;
    }
};

template<RegionTypes RT>
void Region<RT>::info() {
    logging::info(true, "Region: ", this->name);
    logging::info(true, "   Type: ", RT);
    logging::info(true, "   Size: ", this->size());
    logging::info(true, "   IDs : ");
    for (size_t i = 0; i < std::min((size_t) 4, this->size()); i++) {
        logging::info(true, "      ", this->at(i));
    }
}

using NodeRegion    = Region<RegionTypes::NODE>;
using ElementRegion = Region<RegionTypes::ELEMENT>;
using SurfaceRegion = Region<RegionTypes::SURFACE>;

}

#endif    // FEMASTER_REGION_H
