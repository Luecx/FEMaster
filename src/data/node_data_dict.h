//
// Created by f_eggers on 08.11.2024.
//

#ifndef NODE_DATA_DICT_H
#define NODE_DATA_DICT_H

#include "data_storage.h"
#include "dict.h"
#include "../core/core.h"

namespace fem::model {

enum NodeDataEntries : Index {
    POSITION,
    TEMPERATURE,
};

using NodeDataDict  = DataStorage<NodeData>;
using NodeFieldDict = Dict<Dict<NodeData>>;

}



#endif //NODE_DATA_DICT_H
