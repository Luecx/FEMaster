//
// Created by f_eggers on 08.11.2024.
//

#ifndef ELEM_DATA_DICT_H
#define ELEM_DATA_DICT_H


#include "data_storage.h"
#include "dict.h"

#include "../core/types_eig.h"

namespace fem::model {
enum ElementDataEntries : Index {
    TOPO_STIFFNESS,
    TOPO_ANGLES,
};

struct ElementField : NodeData, Namable {
    ElementField(const std::string& name) : Namable(name) {}
    using Ptr = std::shared_ptr<ElementField>;
};


using ElemDataDict  = DataStorage<ElementData>;
using ElemFieldDict = Dict<Dict<ElementField>, ElementDataEntries>;

}



#endif //ELEM_DATA_DICT_H
