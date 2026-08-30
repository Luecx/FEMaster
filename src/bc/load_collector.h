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
 */
struct LoadCollector : model::Collection<Neumann::Ptr> {
    using Ptr = std::shared_ptr<LoadCollector>;
    using model::Collection<Neumann::Ptr>::add;

    explicit LoadCollector(const std::string& name);
    ~LoadCollector() = default;

    /**
     * @brief Applies all conditions in insertion order.
     *
     * Existing structural callers may omit `system_dof_ids` and `lhs`; all
     * purely mechanical loads operate on the RHS only. Thermal analyses that may
     * contain convection pass both pointers so matrix and RHS contributions are
     * assembled through the same condition interface.
     *
     * @param model_data Compiled model data shared by all contained conditions.
     * @param rhs Nodal right-hand-side field receiving contributions.
     * @param time Analysis time used for amplitude evaluation.
     * @param system_dof_ids Optional global system DOF map.
     * @param lhs Optional sparse triplet list receiving LHS contributions.
     */
    void apply(model::ModelData&       model_data,
               model::Field&           rhs,
               Precision               time,
               const SystemDofIds*      system_dof_ids = nullptr,
               TripletList*             lhs = nullptr);

    const std::vector<Neumann::Ptr>& entries() const { return this->_data; }
};

} // namespace fem::bc
