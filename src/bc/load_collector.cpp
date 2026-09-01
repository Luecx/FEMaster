/**
 * @file load_collector.cpp
 * @brief Implements load collector construction and generic RHS dispatch.
 *
 * The collector preserves insertion order and applies every non-null load through
 * the common RHS-only interface. This implementation remains intentionally
 * simple: physical assembly belongs to each concrete boundary condition, while
 * thermal type discrimination between Neumann and Robin conditions is performed
 * by the model's dedicated thermal load path.
 *
 * @see load_collector.h
 * @see model::Model::build_thermal_load_matrix
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "load_collector.h"

namespace fem::bc {

/**
 * Constructs an empty named load collection.
 *
 * @param name Collector name used by loadcase selection and model diagnostics.
 */
LoadCollector::LoadCollector(const std::string& name)
    : model::Collection<Load::Ptr>(name) {}

/**
 * Applies all compatible loads in the collector to a common nodal RHS field.
 *
 * Entries are evaluated in insertion order at the supplied analysis time. Null
 * entries are ignored. The method deliberately uses the generic `Load::apply()`
 * overload and therefore represents an RHS-only dispatch path; a Robin entry
 * fails through its finalized base overload instead of silently omitting its
 * operator contribution.
 *
 * @param model_data Compiled model data required by concrete loads.
 * @param rhs Nodal right-hand-side field receiving accumulated contributions.
 * @param time Analysis time passed to load amplitudes.
 */
void LoadCollector::apply(model::ModelData& model_data,
                          model::Field&     rhs,
                          Precision         time) {
    // Preserve collector order while delegating the physical load definition,
    // amplitude evaluation and geometric integration to each concrete condition.
    for (const auto& load : this->_data) {
        if (load) load->apply(model_data, rhs, time, false);
    }
}

} // namespace fem::bc
