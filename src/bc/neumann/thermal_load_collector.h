/**
 * @file thermal_load_collector.h
 * @brief Declares named collections of thermal Neumann conditions.
 */

#pragma once

#include "thermal_load.h"
#include "../../data/collection.h"

#include <string>
#include <vector>

namespace fem::bc {

struct ThermalLoadCollector : model::Collection<ThermalLoad::Ptr> {
    using Ptr = std::shared_ptr<ThermalLoadCollector>;
    using model::Collection<ThermalLoad::Ptr>::add;

    explicit ThermalLoadCollector(const std::string& name)
        : model::Collection<ThermalLoad::Ptr>(name) {}

    const std::vector<ThermalLoad::Ptr>& entries() const { return this->_data; }
};

} // namespace fem::bc
