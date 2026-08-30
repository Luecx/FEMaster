/**
 * @file load_collector.h
 * @brief Declares named collections of polymorphic load-side boundary conditions.
 *
 * A load collector groups boundary-condition definitions selected together by a
 * loadcase. The common storage type is `Load::Ptr`, allowing structural Neumann,
 * thermal Neumann and Robin conditions to share one named collection while
 * preserving their concrete runtime type for analysis-specific dispatch.
 *
 * The generic `apply()` function is intentionally an RHS-only path retained for
 * structural load assembly. Robin conditions explicitly reject that interface,
 * preventing a mixed boundary condition from losing its operator contribution.
 * Thermal analyses therefore inspect `entries()` and dispatch thermal Neumann and
 * Robin conditions separately through the model-level thermal assembly routine.
 *
 * @see Load
 * @see Neumann
 * @see Robin
 * @see Model::build_thermal_load_matrix
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "load.h"
#include "../data/collection.h"

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Named heterogeneous collection of load-side boundary conditions.
 *
 * The collector inherits the common named `model::Collection` storage and owns
 * shared pointers to load definitions in insertion order. Its generic application
 * path invokes the common RHS-only `Load::apply()` interface and is therefore
 * suitable for structural Neumann collectors. Thermal assembly uses the same
 * stored entries but performs explicit type dispatch because Robin conditions
 * additionally require a system operator and scalar DOF mapping.
 */
struct LoadCollector : model::Collection<Load::Ptr> {
    // Shared ownership type used by ModelData's named collector registry.
    using Ptr = std::shared_ptr<LoadCollector>;

    // Expose the base collection insertion operation for heterogeneous `Load`
    // pointers without adding another forwarding wrapper.
    using model::Collection<Load::Ptr>::add;

    // Construction creates an empty named load collection. The name is retained
    // by the common collection base for selection and diagnostics.
    explicit LoadCollector(const std::string& name);
    ~LoadCollector() = default;

    // Apply every non-null entry through the common RHS-only load interface. This
    // path is used by structural load assembly; Robin entries fail explicitly so
    // their matrix contribution cannot be omitted silently.
    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time);

    // Provide read-only ordered access for analysis-specific dispatch, notably the
    // thermal path that distinguishes `ThermalNeumann` from `Robin` conditions.
    const std::vector<Load::Ptr>& entries() const { return this->_data; }
};

} // namespace fem::bc
