/******************************************************************************
 * @file elem_data_dict.h
 * @brief Declares containers that manage element-level data fields.
 *
 * Provides the `ElementField` helper and storage aliases that mirror the nodal
 * data containers for element-specific matrices such as stiffness updates or
 * orientation angles.
 *
 * @see src/data/elem_data_dict.cpp
 * @see src/data/data_storage.h
 * @see src/data/dict.h
 ******************************************************************************/

#pragma once

#include "../core/types_eig.h"
#include "data_storage.h"
#include "dict.h"
#include "namable.h"

#include <memory>
#include <string>

namespace fem {
namespace model {

/// Enumerates default element data entries consumed by solvers and readers.
enum ElementDataEntries : Index {
    TOPO_STIFFNESS, ///< Scalar density or stiffness scaling factors.
    TOPO_ANGLES     ///< Orientation angles for topology optimisation workflows.
};

/******************************************************************************
 * @struct ElementField
 * @brief Extends `NodeData` with a stable name for element-level datasets.
 ******************************************************************************/
struct ElementField : NodeData, Namable {
    using Ptr = std::shared_ptr<ElementField>; ///< Shared pointer alias for element fields.

    /// Constructs the field with the provided name.
    explicit ElementField(const std::string& name)
        : Namable(name) {}
};

using ElemDataDict = DataStorage<ElementData>;                         ///< Storage for element matrices.
using ElemFieldDict = Dict<Dict<ElementField>, ElementDataEntries>;    ///< Nested dictionary for named element fields.

} // namespace model
} // namespace fem

