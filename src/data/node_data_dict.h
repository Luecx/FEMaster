//
// Created by f_eggers on 08.11.2024.
//

#ifndef NODE_DATA_DICT_H
#define NODE_DATA_DICT_H

#include "data_storage.h"
#include "dict.h"
#include "../core/core.h"

namespace fem::model {

struct NodeField : NodeData {
    using Ptr = std::shared_ptr<NodeField>;

    NodeField() = default;
    NodeField(const NodeData& data) : NodeData(data) {}
    NodeField(const NodeField& field) : NodeData(field) {}
    NodeField(NodeField&& field) noexcept : NodeData(std::move(field)) {}
    NodeField& operator=(const NodeField& field) {
        NodeData::operator=(field);
        return *this;
    }
    NodeField(Index rows, Index cols) : NodeData(rows, cols) {}
};

enum NodeDataEntries : Index {
    POSITION,
    TEMPERATURE,
};

using NodeDataDict  = DataStorage<NodeData>;
using NodeFieldDict = Dict<Dict<NodeField>, NodeDataEntries>;

}



#endif //NODE_DATA_DICT_H
