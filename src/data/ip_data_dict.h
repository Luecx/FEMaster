/******************************************************************************
 * @file ip_data_dict.h
 * @brief Declares integration-point data containers.
 *
 * Integration-point fields reuse the generic `DataStorage` infrastructure to
 * keep track of stress, strain, or other state variables evaluated at Gauss
 * points.
 *
 * @see src/data/ip_data_dict.cpp
 * @see src/data/data_storage.h
 ******************************************************************************/

#pragma once

#include "../core/types_eig.h"
#include "data_storage.h"
#include "dict.h"

namespace fem {
namespace model {

/// Enumerates standard integration-point data entries.
enum IPDataEntries : Index {
    STRESS ///< Cauchy stress tensor stored at each integration point.
};

using IPDataDict = DataStorage<IPData>; ///< Storage alias for integration-point state variables.

} // namespace model
} // namespace fem

