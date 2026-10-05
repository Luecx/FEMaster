/**
 * @file region.h
 * @brief Defines typed named collections of model entity identifiers.
 *
 * The model-data subsystem uses Region to group node, element, surface and line
 * identifiers under an immutable name. Collection supplies storage, insertion
 * policies and parent forwarding; the compile-time region kind supplies semantic
 * identity and diagnostic output. Sets manages named registrations and aggregates.
 *
 * Identifier interpretation and remapping belong to the model and instance
 * compilation. Region does not resolve or own the referenced entities.
 *
 * @see Region
 * @see Collection
 * @see Sets
 * @see RegionTypes
 *
 * @author Finn Eggers
 * @date 05.10.2026
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
 * @brief Names and categorizes a collection of model entity identifiers.
 *
 * The RegionTypes template argument fixes the entity kind without adding runtime
 * storage. The inherited Collection<ID> owns the identifier values and immutable
 * name. Construction preserves insertion order and allows repeated identifiers;
 * consumers may change those inherited policies where their use requires it.
 *
 * Entity ownership, identifier validity and part-to-assembly remapping remain
 * model responsibilities. Parent links and shared ownership follow Collection's
 * contract. Diagnostics report the kind and size; info() prints at most four IDs,
 * while the member stream operation emits metadata and an IDs label only.
 *
 * @tparam RT NODE, ELEMENT, SURFACE or LINE entity kind.
 */
template<RegionTypes RT>
struct Region : public Collection<ID> {
    // Shared ownership of this concrete entity-kind collection.
    using Ptr = std::shared_ptr<Region<RT>>;

    // Construction with insertion-order storage and permitted duplicates.
    explicit Region(std::string name)
        : Collection<ID>(std::move(name), true, false) {}

    // Log region metadata and up to four stored identifiers.
    void info();

    // Append metadata and the IDs label to the stream; identifier values are omitted.
    std::ostream& operator<<(std::ostream& os) const {
        os << "Region: " << this->name;
        os << "   Type: " << RT;
        os << "   Size: " << this->size();
        os << "   IDs : ";
        return os;
    }
};

/**
 * @brief Logs the region kind, size and a short identifier preview.
 *
 * Print the immutable name and numeric compile-time region kind, then report the
 * stored size and at most the first four identifiers in the current storage
 * order. This operation neither sorts nor changes the region or its parent.
 */
template<RegionTypes RT>
void Region<RT>::info() {
    // Report metadata before limiting the identifier preview to four entries.
    logging::info(true, "Region: ", this->name);
    logging::info(true, "   Type: ", RT);
    logging::info(true, "   Size: ", this->size());
    logging::info(true, "   IDs : ");
    for (size_t i = 0; i < std::min<size_t>(4, this->size()); ++i) {
        logging::info(true, "      ", this->at(i));
    }
}

// Entity-kind aliases retain the same identifier storage and insertion policies.
using NodeRegion = Region<RegionTypes::NODE>;
using ElementRegion = Region<RegionTypes::ELEMENT>;
using SurfaceRegion = Region<RegionTypes::SURFACE>;
using LineRegion = Region<RegionTypes::LINE>;
} // namespace model
} // namespace fem
