/**
 * @file region.h
 * @brief Declares generic region collections used by FEM models.
 *
 * A region groups identifiers (nodes, elements, surfaces) under a common name
 * and provides helper functions for logging and iteration.
 *
 * @see src/data/region.cpp
 * @see src/data/region_type.h
 * @see src/data/collection.h
 */

#pragma once

#include "collection.h"
#include "region_type.h"

#include "../core/logging.h"
#include "../core/types_num.h"

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>

namespace fem {
namespace model {

/**
 * @struct Region
 * @brief Named collection of entity identifiers.
 *
 * @tparam RT Compile-time region type (`NODE`, `ELEMENT`, or `SURFACE`).
 */
template<RegionTypes RT>
struct Region : public Collection<ID> {
    using Ptr = std::shared_ptr<Region<RT>>; ///< Shared pointer alias for region ownership.

    /**
     * @brief Constructs a region with duplicate tracking but without sorting.
     */
    explicit Region(std::string name)
        : Collection<ID>(std::move(name), true, false) {}

    /**
     * @brief Emits logging information about the region contents.
     */
    void info();

    /**
     * @brief Streams a textual representation into `os`.
     */
    std::ostream& operator<<(std::ostream& os) const {
        os << "Region: " << this->name;
        os << "   Type: " << RT;
        os << "   Size: " << this->size();
        os << "   IDs : ";
        return os;
    }
};

/**
 * @copydoc Region<RT>::info
 */
template<RegionTypes RT>
void Region<RT>::info() {
    logging::info(true, "Region: ", this->name);
    logging::info(true, "   Type: ", RT);
    logging::info(true, "   Size: ", this->size());
    logging::info(true, "   IDs : ");
    for (size_t i = 0; i < std::min<size_t>(4, this->size()); ++i) {
        logging::info(true, "      ", this->at(i));
    }
}

using NodeRegion = Region<RegionTypes::NODE>;       ///< Region alias for nodes.
using ElementRegion = Region<RegionTypes::ELEMENT>; ///< Region alias for elements.
using SurfaceRegion = Region<RegionTypes::SURFACE>; ///< Region alias for surfaces.
using LineRegion = Region<RegionTypes::LINE>;       ///< Region alias for lines (1D geometry).

} // namespace model
} // namespace fem
