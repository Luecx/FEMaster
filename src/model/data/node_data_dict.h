//
// Created by f_eggers on 08.11.2024.
//

#ifndef NODE_DATA_DICT_H
#define NODE_DATA_DICT_H

#include "../collection/data_storage.h"
#include "../../core/core.h"

namespace fem::model {

enum NodalDataEntries : Index {
    POSITION
};

using NodeDataDict = DataStorage<NodeData>;

}



#endif //NODE_DATA_DICT_H
