/**
 * @file node_data_dict.h
 * @brief Declares containers that store nodal field data.
 *
 * Provides the `NodeField` helper as well as type aliases that combine
 * `DataStorage` and `Dict` to manage named nodal matrices.
 *
 * @see src/data/node_data_dict.cpp
 * @see src/data/data_storage.h
 * @see src/data/dict.h
 */

#pragma once

#include "../core/types_eig.h"
#include "data_storage.h"
#include "dict.h"
#include "namable.h"

#include <memory>
#include <string>
#include <utility>

namespace fem {
namespace model {

/**
 * @struct NodeField
 * @brief Wraps nodal data together with a human-readable name.
 */
struct NodeField : Namable, NodeData {
    using Ptr = std::shared_ptr<NodeField>; ///< Shared pointer alias for node fields.

    /// Constructs an empty field with the provided name.
    explicit NodeField(const std::string& name)
        : Namable(name) {}

    /// Copy constructor required by `Dict` when cloning entries.
    NodeField(const NodeField& field)
        : Namable(field), NodeData(field) {}

    /// Move constructor used when transferring ownership between containers.
    NodeField(NodeField&& field) noexcept
        : Namable(std::move(field)), NodeData(std::move(field)) {}

    /// Copy assignment to replace the underlying node data while keeping the name.
    NodeField& operator=(const NodeField& field) {
        NodeData::operator=(field);
        return *this;
    }

    /// Constructs a nameless field with custom dimensions.
    NodeField(Index rows, Index cols)
        : Namable(""), NodeData(rows, cols) {}

    /// Constructs a named field with custom dimensions.
    NodeField(const std::string& name, Index rows, Index cols)
        : Namable(name), NodeData(rows, cols) {}
};

/// Enumerates standard nodal data entries used throughout the solver.
enum NodeDataEntries : Index {
    POSITION,   ///< Cartesian coordinates of each node.
    TEMPERATURE ///< Temperature degree of freedom per node.
};

using NodeDataDict = DataStorage<NodeData>;                           ///< Storage for nodal matrices.
using NodeFieldDict = Dict<Dict<NodeField>, NodeDataEntries>;          ///< Nested dictionary for named nodal fields.

} // namespace model
} // namespace fem

