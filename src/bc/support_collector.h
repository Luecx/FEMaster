/**
 * @file support_collector.h
 * @brief Declares named heterogeneous support collections.
 *
 * A support collector may own mechanical displacement constraints and thermal
 * temperature constraints. Each load case collects only its compatible concrete
 * support type so local degree-of-freedom identifiers are never interpreted by
 * the wrong solver system.
 *
 * @see SupportInterface
 * @see Support
 * @see Temperature
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../data/collection.h"
#include "support_interface.h"

#include <type_traits>

namespace fem {
namespace model {
struct ModelData;
}

namespace bc {

/**
 * @brief Owns prescribed boundary conditions under one shared collector name.
 *
 * Mechanical and thermal support definitions coexist as polymorphic objects in
 * insertion order. A requesting load case selects one concrete support type when
 * it asks for equations; other definitions remain in the collector but do not
 * contribute. This permits one input-level support namespace without confusing
 * structural degree of freedom zero with the scalar thermal temperature degree
 * of freedom.
 */
struct SupportCollector : model::Collection<SupportInterface::Ptr> {
    using Ptr = std::shared_ptr<SupportCollector>;
    using model::Collection<SupportInterface::Ptr>::add;

    explicit SupportCollector(const std::string& name);
    ~SupportCollector() = default;

    // Apply entries of the requested concrete support type in insertion order.
    // Mechanical and thermal load cases use this filter with their respective
    // support implementation so a shared collector remains solver-safe.
    template<typename SupportType>
    constraint::Equations get_equations(model::ModelData& model_data) {
        // Restrict selection to concrete implementations of the common support
        // contract at compile time.
        static_assert(std::is_base_of_v<SupportInterface, SupportType>,
            "SupportType must derive from SupportInterface");

        // Preserve collector order while applying only the requested support
        // implementation.
        constraint::Equations equations{};
        for (const auto& support : this->_data) {
            if (auto typed_support = std::dynamic_pointer_cast<SupportType>(support)) {
                typed_support->apply(model_data, equations);
            }
        }
        return equations;
    }

    const std::vector<SupportInterface::Ptr>& entries() const { return this->_data; }
};

} // namespace bc
} // namespace fem
