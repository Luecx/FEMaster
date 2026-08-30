/**
 * @file load_collector.cpp
 * @brief Implements collective application of load-side boundary conditions.
 *
 * Each contained condition receives the same right-hand-side field and optional
 * system assembly objects. This keeps load collection independent of the exact
 * algebraic contribution made by an individual condition: conventional loads
 * modify only RHS, while convection may additionally append LHS triplets.
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "load_collector.h"

#include "../data/field.h"

namespace fem::bc {

LoadCollector::LoadCollector(const std::string& name)
    : model::Collection<Neumann::Ptr>(name) {}

/**
 * Applies every stored condition in insertion order.
 *
 * @param model_data Compiled model data shared by all conditions.
 * @param rhs Nodal right-hand-side field receiving load contributions.
 * @param time Analysis time used for amplitude evaluation.
 * @param system_dof_ids Optional global system DOF map for matrix-producing loads.
 * @param lhs Optional sparse triplet list for matrix-producing loads.
 */
void LoadCollector::apply(model::ModelData&       model_data,
                          model::Field&           rhs,
                          Precision               time,
                          const SystemDofIds*      system_dof_ids,
                          TripletList*             lhs) {
    for (const auto& load : this->_data) {
        if (load) {
            load->apply(model_data, rhs, time, false, system_dof_ids, lhs);
        }
    }
}

} // namespace fem::bc
