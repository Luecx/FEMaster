//
// Created by f_eggers on 08.11.2024.
//

#ifndef NODE_DATA_DICT_H
#define NODE_DATA_DICT_H

#include "../core/core.h"
#include "data_storage.h"
#include "dict.h"
#include "namable.h"

namespace fem::model {

struct NodeField : Namable, NodeData{
    using Ptr = std::shared_ptr<NodeField>;

    NodeField(const std::string& name) : Namable(name) {};
    NodeField(const NodeField& field) : Namable(field), NodeData(field) {}
    NodeField(NodeField&& field) noexcept : Namable(std::move(field)), NodeData(std::move(field)) {}
    NodeField& operator=(const NodeField& field) {
        NodeData::operator=(field);
        return *this;
    }
    NodeField(Index rows, Index cols) : Namable(""), NodeData(rows, cols) {}
    NodeField(const std::string& name, Index rows, Index cols) : Namable(name), NodeData(rows, cols) {}
};

enum NodeDataEntries : Index {
    POSITION,
    TEMPERATURE,
};

using NodeDataDict  = DataStorage<NodeData>;
using NodeFieldDict = Dict<Dict<NodeField>, NodeDataEntries>;

}



#endif //NODE_DATA_DICT_H
