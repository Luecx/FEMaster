/**
 * @file load_collector.cpp
 * @brief Implements superposed RHS and mixed-operator assembly for load collectors.
 *
 * The collector itself contains no finite-element integration formulas. It
 * traverses the stored polymorphic loads and delegates each physical
 * contribution to the concrete condition. RHS assembly is applied to every
 * entry, while matrix assembly selects only conditions from the `Mixed`
 * category.
 *
 * This separation mirrors the discrete equilibrium structure: Neumann and mixed
 * conditions may both contribute to the external generalized load vector, but
 * only mixed conditions provide an additional operator term.
 *
 * @see LoadCollector
 * @see Load
 * @see Mixed
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#include "load_collector.h"
#include "mixed/mixed.h"

#include "../data/field.h"

namespace fem::bc {

LoadCollector::LoadCollector(const std::string& name)
    : model::Collection<Load::Ptr>(name) {}

/**
 * Superimposes all stored right-hand-side contributions.
 *
 * Every load receives the same global nodal field and adds its equivalent
 * generalized force contribution in place. If the individual loads are denoted
 * by `f_a`, the resulting operation is
 *
 *     rhs <- rhs + sum_a f_a.
 *
 * No temporary collector-level field is required because finite-element load
 * assembly is additive by construction. Coordinate transformations, geometric
 * integration and amplitude evaluation remain local to the concrete loads.
 *
 * @param model_data Model data required by the concrete load implementations.
 * @param rhs Generalized nodal right-hand-side field modified in place.
 * @param time Analysis time forwarded to amplitude-dependent conditions.
 */
void LoadCollector::apply(model::ModelData& model_data, model::Field& rhs, Precision time) {
    // Apply every valid load directly to the common RHS so overlapping
    // definitions superimpose naturally at shared nodal degrees of freedom
    for (const auto& load : this->_data) {
        if (!load) {
            continue;
        }
        load->apply(model_data, rhs, time);
    }
}

/**
 * Assembles the operator terms of mixed conditions stored in the collector.
 *
 * Pure Neumann entries are right-hand-side terms and therefore do not
 * participate in this pass. A mixed condition contributes sparse triplets in
 * active system numbering, schematically
 *
 *     K_b <- K_b + K_a.
 *
 * The caller owns the final sign and placement of the accumulated boundary
 * operator in the global residual or tangent matrix. This function only
 * superimposes the triplets supplied by each mixed condition.
 *
 * @param model_data Model data required by the mixed-condition implementation.
 * @param system_dof_ids Mapping from model DOFs to active system equation IDs.
 * @param matrix Sparse triplet list receiving the operator contributions.
 * @param time Analysis time forwarded to amplitude-dependent conditions.
 * @param ignore_amplitude Evaluate nominal condition magnitudes with unit
 *                         amplitude when true.
 */
void LoadCollector::apply_matrix(model::ModelData&   model_data,
                                 const SystemDofIds& system_dof_ids,
                                 TripletList&        matrix,
                                 Precision           time,
                                 bool                ignore_amplitude) {
    // Select only Mixed entries from the shared Load collection. The dynamic
    // cast expresses the algebraic category explicitly without requiring the
    // collector to know any concrete Robin or convection type.
    for (const auto& load : this->_data) {
        if (!load) {
            continue;
        }

        auto mixed = std::dynamic_pointer_cast<Mixed>(load);
        if (!mixed) {
            continue;
        }

        // Let the concrete mixed condition map its model-level contribution to
        // the active sparse system and append the resulting triplets
        mixed->apply_matrix(model_data, system_dof_ids, matrix, time, ignore_amplitude);
    }
}

} // namespace fem::bc
