/**
 * @file load_collector.h
 * @brief Declares a named collection of load-side boundary conditions.
 *
 * Structural loads, prescribed thermal heat flux and convection share the
 * `Neumann` abstraction and retain insertion order inside one named collector.
 * The collector forwards one common assembly context to every contained load.
 * Pure RHS conditions ignore the optional system objects, whereas conditions
 * such as convection use them to add LHS matrix entries in the same pass.
 *
 * @see Neumann
 * @see HeatFlux
 * @see Convection
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
 * @brief Owns load-side boundary conditions under one shared collector name.
 *
 * Entries retain input order and are also propagated to the aggregate load
 * collector maintained by `model::Sets`. `apply()` forwards one common model,
 * time and algebraic assembly context to every condition. This allows pure
 * Neumann loads to contribute only to the nodal RHS while Robin conditions can
 * append matrix terms without a separate collector abstraction.
 */
struct LoadCollector : model::Collection<Neumann::Ptr> {
    using Ptr = std::shared_ptr<LoadCollector>;
    using model::Collection<Neumann::Ptr>::add;

    explicit LoadCollector(const std::string& name);
    ~LoadCollector() = default;

    // Apply all conditions in insertion order. Structural callers may omit the
    // optional system objects; thermal convection requires both for its LHS.
    void apply(model::ModelData&       model_data,
               model::Field&           rhs,
               Precision               time,
               const SystemDofIds*      system_dof_ids = nullptr,
               TripletList*             lhs = nullptr);

    const std::vector<Neumann::Ptr>& entries() const { return this->_data; }
};

} // namespace fem::bc
