/**
 * @file load_collector.h
 * @brief Declares named collections of load-side boundary conditions.
 */

#pragma once

#include "load.h"
#include "../data/collection.h"

namespace fem::bc {

/**
 * @brief Owns heterogeneous load-side conditions while preserving insertion order.
 *
 * Model-level assembly selects the required derived interface explicitly:
 * structural/thermal Neumann conditions contribute only to RHS fields, while
 * Robin conditions additionally emit symbolic equation rows.
 */
struct LoadCollector : model::Collection<Load::Ptr> {
    using Ptr = std::shared_ptr<LoadCollector>;
    using model::Collection<Load::Ptr>::add;

    explicit LoadCollector(const std::string& name);
    ~LoadCollector() = default;

    const std::vector<Load::Ptr>& entries() const { return this->_data; }
};

} // namespace fem::bc
