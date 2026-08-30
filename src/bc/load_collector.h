/**
 * @file load_collector.h
 * @brief Declares a named collection of natural boundary conditions.
 *
 * Structural loads and thermal natural boundary conditions share the `Neumann`
 * abstraction and retain insertion order inside one named collector.
 *
 * @see Neumann
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "neumann/neumann.h"
#include "../data/collection.h"

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Owns Neumann conditions under one shared load-collector name.
 */
struct LoadCollector : model::Collection<Neumann::Ptr> {
    using Ptr = std::shared_ptr<LoadCollector>;
    using model::Collection<Neumann::Ptr>::add;

    explicit LoadCollector(const std::string& name);
    ~LoadCollector() = default;

    // Apply all conditions in insertion order to the supplied nodal field.
    void apply(model::ModelData& model_data, model::Field& bc, Precision time);

    const std::vector<Neumann::Ptr>& entries() const { return this->_data; }
};

} // namespace fem::bc
