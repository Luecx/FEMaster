/**
 * @file load_collector.h
 * @brief Declares a named collection of polymorphic loads.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../data/collection.h"
#include "load.h"

namespace fem {
namespace model {
class ModelData;
}

namespace bc {

struct LoadCollector : model::Collection<Load::Ptr> {
    using Ptr = std::shared_ptr<LoadCollector>;
    using model::Collection<Load::Ptr>::add;

    explicit LoadCollector(const std::string& name);
    ~LoadCollector() = default;

    void apply(model::ModelData& model_data, model::Field& bc, Precision time);

    const std::vector<Load::Ptr>& entries() const { return this->_data; }
};

} // namespace bc
} // namespace fem
