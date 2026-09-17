/**
 * @file load_collector.h
 * @brief Defines named collections of right-hand-side-capable conditions.
 *
 * A load collector groups Neumann and mixed conditions that belong to the same
 * analysis load definition. The collector superimposes their equivalent nodal
 * right-hand-side contributions and, when requested, extracts the additional
 * operator terms provided by mixed conditions.
 *
 * Keeping the collector on the shared `Load` capability allows the load-case
 * layer to assemble one named collection without depending on the concrete
 * mathematical category of every entry.
 *
 * @see Load
 * @see Neumann
 * @see Mixed
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include "load.h"
#include "../data/collection.h"

#include <memory>
#include <string>
#include <vector>

namespace fem::bc {

/**
 * @brief Named superposition of load-like boundary conditions.
 *
 * If a collector contains conditions `a = 1, ..., n`, RHS assembly performs the
 * linear accumulation
 *
 *     f_ext <- f_ext + sum_a f_a.
 *
 * Mixed conditions may additionally contribute sparse operator terms
 *
 *     K_b <- K_b + sum_a K_a,
 *
 * while pure Neumann conditions are ignored by matrix assembly. The collector
 * does not own any finite-element integration logic itself; every condition
 * computes its own physically consistent nodal or matrix contribution.
 */
struct LoadCollector : model::Collection<Load::Ptr> {
    // Types
    using Ptr = std::shared_ptr<LoadCollector>;
    using model::Collection<Load::Ptr>::add;

    // Construction
    explicit LoadCollector(const std::string& name);
    ~LoadCollector() = default;

    // Superimpose the RHS contribution of every stored load into the supplied
    // generalized nodal field.
    void apply(model::ModelData& model_data, model::Field& rhs, Precision time);

    // Assemble only the operator contributions of stored Mixed conditions in
    // active system numbering. Pure Neumann entries contribute no matrix terms.
    void apply_matrix(model::ModelData&   model_data,
                      const SystemDofIds& system_dof_ids,
                      TripletList&        matrix,
                      Precision           time,
                      bool                ignore_amplitude = false);

    // Read-only access to the stored polymorphic load definitions
    const std::vector<Load::Ptr>& entries() const { return this->_data; }
};

} // namespace fem::bc
