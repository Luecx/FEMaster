//
// Created by f_eggers on 08.11.2024.
//

#ifndef ELEM_DATA_DICT_H
#define ELEM_DATA_DICT_H


#include "../collection/data_storage.h"
#include "../../core/core.h"

namespace fem::model {
enum ElementDataEntries : Index {
    TOPO_DENSITY,
    MAT_ANGLES,
};

using ElemDataDict = DataStorage<ElementData>;

}



#endif //ELEM_DATA_DICT_H
