/**
 * @file thermal_collector.cpp
 * @brief Implements algebraically separated assembly of thermal BC collections.
 *
 * Thermal conditions are stored by physical domain but assembled by mathematical
 * role. Dirichlet entries generate scalar temperature constraints, Neumann and
 * Mixed entries contribute to the nodal heat-flow RHS, and only Mixed entries
 * contribute an additional operator.
 *
 * This keeps the eventual thermal loadcase simple: it can request all three
 * contributions from one named collector without knowing concrete condition
 * types, while the order of system construction remains
 *
 *     physical operator and RHS -> Dirichlet constraint transformation.
 *
 * @see ThermalCollector
 * @see ThermalCondition
 *
 * @author Finn Eggers
 * @date 18.09.2026
 */

#include "thermal_collector.h"

#include "dirichlet/dirichlet.h"
#include "load.h"
#include "mixed/mixed.h"

#include "../core/logging.h"

namespace fem::bc {

ThermalCollector::ThermalCollector(const std::string& name)
    : model::Collection<ThermalCondition::Ptr>(name) {}

/**
 * Adds one thermal boundary condition to the collection.
 *
 * The thermal marker identifies the physical field but does not itself define an
 * algebraic contribution. Every stored object must therefore also implement
 * either the Dirichlet interface or the Load interface. Mixed conditions satisfy
 * the latter through their inherited RHS contract and additionally expose their
 * operator contribution during `apply_matrix()`.
 *
 * @param condition Thermal boundary condition to store.
 */
void ThermalCollector::add(const ThermalCondition::Ptr& condition) {
    // Reject null definitions before they can enter the named collector
    logging::error(condition != nullptr,
        "THERMAL: cannot add a null boundary condition");

    // Verify that the physical-domain marker is paired with a supported
    // algebraic category. A thermal condition is either prescribed through
    // C T = d or contributes to the assembled thermal balance.
    const bool is_dirichlet = std::dynamic_pointer_cast<Dirichlet>(condition) != nullptr;
    const bool is_load      = std::dynamic_pointer_cast<Load>(condition)      != nullptr;

    logging::error(is_dirichlet || is_load,
        "THERMAL: condition has no Dirichlet, Neumann or Mixed assembly interface");

    // Preserve polymorphic ownership through the thermal marker
    model::Collection<ThermalCondition::Ptr>::add(condition);
}

/**
 * Collects the Dirichlet equations of all prescribed thermal primary variables.
 *
 * Only thermal objects that also derive from `Dirichlet` participate. For the
 * current condition set this means prescribed temperature rows
 *
 *     T_i = T_bar.
 *
 * Heat-flux and convection entries do not generate algebraic constraints and are
 * ignored by this pass.
 *
 * @param model_data Compiled model data used by concrete Dirichlet conditions.
 * @return Concatenated thermal constraint equations.
 */
constraint::Equations ThermalCollector::get_equations(model::ModelData& model_data) {
    constraint::Equations equations{};

    // Append only essential thermal boundary conditions to C T = d
    for (const auto& condition : this->_data) {
        if (!condition) {
            continue;
        }

        auto dirichlet = std::dynamic_pointer_cast<Dirichlet>(condition);
        if (!dirichlet) {
            continue;
        }

        dirichlet->apply(model_data, equations);
    }

    return equations;
}

/**
 * Superimposes every load-like thermal contribution into the scalar nodal RHS.
 *
 * Neumann conditions contribute prescribed heat input directly. Mixed conditions
 * contribute only their prescribed source part during this pass; their
 * temperature-dependent operator is assembled separately by `apply_matrix()`.
 * Consequently the accumulated field represents
 *
 *     q_total = sum q_N + sum q_M.
 *
 * @param model_data Compiled model data required by concrete conditions.
 * @param rhs Scalar nodal thermal right-hand-side field modified in place.
 * @param time Analysis time forwarded to optional amplitudes.
 * @param ignore_amplitude Use nominal condition values when true.
 */
void ThermalCollector::apply_rhs(model::ModelData& model_data,
                                 model::Field&     rhs,
                                 Precision         time,
                                 bool              ignore_amplitude) {
    // Load is the shared RHS capability of Neumann and Mixed conditions
    for (const auto& condition : this->_data) {
        if (!condition) {
            continue;
        }

        auto load = std::dynamic_pointer_cast<Load>(condition);
        if (!load) {
            continue;
        }

        load->apply(model_data, rhs, time, ignore_amplitude);
    }
}

/**
 * Assembles all unknown-dependent thermal boundary operators.
 *
 * Only Mixed conditions participate. For linear convection this pass adds
 *
 *     K_h = integral_Gamma h N N^T dGamma
 *
 * to the supplied sparse triplet list. Dirichlet conditions remain algebraic
 * constraints and pure Neumann conditions have no operator contribution.
 *
 * @param model_data Compiled model data required by concrete conditions.
 * @param system_dof_ids Scalar node-to-active-temperature equation mapping.
 * @param matrix Sparse triplet list receiving mixed boundary operators.
 * @param time Analysis time forwarded to optional amplitudes.
 * @param ignore_amplitude Use nominal condition values when true.
 */
void ThermalCollector::apply_matrix(model::ModelData&   model_data,
                                    const SystemDofIds& system_dof_ids,
                                    TripletList&        matrix,
                                    Precision           time,
                                    bool                ignore_amplitude) {
    // Select only Mixed thermal conditions for the operator assembly pass
    for (const auto& condition : this->_data) {
        if (!condition) {
            continue;
        }

        auto mixed = std::dynamic_pointer_cast<Mixed>(condition);
        if (!mixed) {
            continue;
        }

        mixed->apply_matrix(model_data, system_dof_ids, matrix, time, ignore_amplitude);
    }
}

} // namespace fem::bc
