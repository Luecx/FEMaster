/**
 * @file support_collector.h
 * @brief Declares named collections of Dirichlet boundary conditions.
 *
 * Support collectors store prescribed primary variables through the common
 * `Dirichlet` interface. Structural and thermal conditions may share a named
 * collector while each analysis selects only the compatible concrete type.
 *
 * @see Dirichlet
 * @see Support
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "dirichlet/dirichlet.h"
#include "../data/collection.h"

#include <type_traits>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Owns prescribed Dirichlet conditions under one collector name.
 *
 * Conditions remain in insertion order. The unqualified `get_equations()` path
 * preserves structural behavior, while typed retrieval lets other analyses
 * select their own Dirichlet condition without another collector hierarchy.
 */
struct SupportCollector : model::Collection<Dirichlet::Ptr> {
    using Ptr = std::shared_ptr<SupportCollector>;
    using model::Collection<Dirichlet::Ptr>::add;

    explicit SupportCollector(const std::string& name);
    ~SupportCollector() = default;

    // Structural compatibility path and typed equation generation
    constraint::Equations get_equations(model::ModelData& model_data);

    template<typename DirichletType>
    constraint::Equations get_equations(model::ModelData& model_data) {
        static_assert(std::is_base_of_v<Dirichlet, DirichletType>,
            "DirichletType must derive from Dirichlet");

        constraint::Equations equations{};
        for (const auto& condition : this->_data) {
            if (auto typed_condition = std::dynamic_pointer_cast<DirichletType>(condition)) {
                typed_condition->apply(model_data, equations);
            }
        }
        return equations;
    }

    const std::vector<Dirichlet::Ptr>& entries() const { return this->_data; }
};

} // namespace fem::bc
