/**
 * @file support_collector.h
 * @brief Declares a named collection of structural supports.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../data/collection.h"
#include "support.h"

namespace fem {
namespace model {
class ModelData;
}

namespace bc {

struct SupportCollector : model::Collection<Support> {
    using Ptr = std::shared_ptr<SupportCollector>;
    using model::Collection<Support>::add;

    explicit SupportCollector(const std::string& name);
    ~SupportCollector() = default;

    constraint::Equations get_equations(model::ModelData& model_data);

    const std::vector<Support>& entries() const { return this->_data; }
};

} // namespace bc
} // namespace fem
