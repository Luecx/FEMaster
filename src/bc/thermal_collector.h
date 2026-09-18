/**
 * @file thermal_collector.h
 * @brief Defines named collections containing all thermal boundary conditions.
 *
 * A thermal analysis needs three algebraically different boundary contributions:
 * prescribed temperatures generate constraint equations, prescribed heat fluxes
 * contribute only to the scalar thermal RHS, and mixed conditions such as
 * convection contribute to both the RHS and thermal operator.
 *
 * `ThermalCollector` stores these definitions together because they belong to
 * one physical thermal boundary-condition set. Assembly remains split into
 * separate passes so the discrete system retains the form
 *
 *     (K_T + K_b) T = q + q_b,
 *
 * followed by enforcement of the collected Dirichlet equations.
 *
 * @see ThermalCondition
 * @see Dirichlet
 * @see Neumann
 * @see Mixed
 *
 * @author Finn Eggers
 * @date 18.09.2026
 */

#pragma once

#include "thermal.h"

#include "../constraints/types/equation.h"
#include "../core/types_eig.h"
#include "../core/types_num.h"
#include "../data/collection.h"
#include "../data/field.h"

#include <memory>
#include <string>
#include <vector>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Named set of thermal Dirichlet, Neumann and mixed conditions.
 *
 * The collector is typed through `ThermalCondition`, so structural supports and
 * structural-only loads cannot be inserted accidentally. Each stored object is
 * still dispatched through its algebraic category:
 *
 * - Dirichlet conditions append rows `C T = d`,
 * - load-like conditions add to the scalar nodal thermal RHS,
 * - Mixed conditions additionally append sparse boundary-operator triplets.
 *
 * The collector performs no surface integration or physical constitutive work
 * itself. Concrete boundary conditions own those formulations.
 */
struct ThermalCollector : model::Collection<ThermalCondition::Ptr> {
    // Types
    using Ptr = std::shared_ptr<ThermalCollector>;

    // Construction
    explicit ThermalCollector(const std::string& name);
    ~ThermalCollector() = default;

    // Insert one thermal condition after validating that it also provides one
    // of the supported algebraic boundary-condition interfaces.
    void add(const ThermalCondition::Ptr& condition);

    // Collect all prescribed-temperature equations C T = d. Flux and mixed
    // conditions do not participate in this pass.
    constraint::Equations get_equations(model::ModelData& model_data);

    // Superimpose all load-like thermal contributions into the scalar nodal RHS.
    // This includes pure heat flux and the prescribed source part of Mixed BCs.
    void apply_rhs(model::ModelData& model_data,
                   model::Field&     rhs,
                   Precision         time,
                   bool              ignore_amplitude = false);

    // Assemble only the unknown-dependent operator terms of Mixed thermal
    // conditions in active scalar temperature-system numbering.
    void apply_matrix(model::ModelData&   model_data,
                      const SystemDofIds& system_dof_ids,
                      TripletList&        matrix,
                      Precision           time,
                      bool                ignore_amplitude = false);

    // Read-only access used by diagnostics and future parser/loadcase reporting
    const std::vector<ThermalCondition::Ptr>& entries() const { return this->_data; }
};

} // namespace fem::bc
