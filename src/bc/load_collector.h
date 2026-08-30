/**
 * @file load_collector.h
 * @brief Declares named collections of load-side boundary conditions.
 */

#pragma once

#include "load.h"
#include "../data/collection.h"

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

struct LoadCollector : model::Collection<Load::Ptr> {
    using Ptr = std::shared_ptr<LoadCollector>;
    using model::Collection<Load::Ptr>::add;

    explicit LoadCollector(const std::string& name);
    ~LoadCollector() = default;

    // Generic RHS-only application retained for structural load paths. Robin
    // entries reject this dispatch and therefore cannot be silently misused.
    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time);

    const std::vector<Load::Ptr>& entries() const { return this->_data; }
};

} // namespace fem::bc
