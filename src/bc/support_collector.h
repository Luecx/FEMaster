/**
 * @file support_collector.h
 * @brief Defines named collections of structural Dirichlet supports.
 *
 * Structural supports are grouped separately from generic loads because they do
 * not assemble generalized forces. Instead each support expands its target
 * region and generates rows of the global constraint system
 *
 *     C u = d.
 *
 * `SupportCollector` concatenates those rows for one named support definition
 * used by an analysis step.
 *
 * @see Support
 * @see Dirichlet
 * @see constraint::Equation
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include "dirichlet/support.h"
#include "../constraints/types/equation.h"
#include "../data/collection.h"

#include <memory>
#include <string>

namespace fem::bc {

/**
 * @brief Named collection of structural support prescriptions.
 *
 * Every stored `Support` may generate one or more scalar equations for each
 * target node. `get_equations()` concatenates these rows without eliminating
 * redundancy; rank analysis and treatment of dependent equations remain the
 * responsibility of the downstream constraint system.
 */
struct SupportCollector : model::Collection<Support> {
    // Types
    using Ptr = std::shared_ptr<SupportCollector>;
    using model::Collection<Support>::add;

    // Construction
    explicit SupportCollector(const std::string& name);
    ~SupportCollector() = default;

    // Expand all stored supports and return the accumulated rows C_i u = d_i
    // in the order in which the supports are stored.
    constraint::Equations get_equations(model::ModelData& model_data);

    // Expose the stored support definitions for read-only diagnostics and model
    // overview output without exposing the protected collection storage itself.
    const std::vector<Support>& entries() const { return this->_data; }
};

} // namespace fem::bc
