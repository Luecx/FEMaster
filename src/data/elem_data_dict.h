//
// Created by f_eggers on 08.11.2024.
//

#ifndef ELEM_DATA_DICT_H
#define ELEM_DATA_DICT_H


#include "data_storage.h"
#include "dict.h"

#include "../core/core.h"

namespace fem::model {
enum ElementDataEntries : Index {
    TOPO_STIFFNESS,
    TOPO_ANGLES,
};

using ElemDataDict  = DataStorage<ElementData>;
using ElemFieldDict = Dict<Dict<ElementData>>;

}



#endif //ELEM_DATA_DICT_H
