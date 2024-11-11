//
// Created by Finn Eggers on 02.11.24.
//

#ifndef COLLECTION_REGION_TYPE_H
#define COLLECTION_REGION_TYPE_H

namespace fem {

using RegionType = int;

enum RegionTypes : RegionType {
    NODE,
    ELEMENT,
    SURFACE,
};

}

#endif //COLLECTION_REGION_TYPE_H
