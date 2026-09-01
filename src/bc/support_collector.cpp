/**
 * @file support_collector.cpp
 * @brief Implements construction and structural dispatch of Dirichlet collectors.
 *
 * The implementation keeps the historical unqualified support path as a thin
 * structural specialization of the typed collector interface. Thermal analyses
 * call the same template with `Temperature`, so both physics share collector
 * ownership while retaining analysis-specific equation selection.
 *
 * @see support_collector.h
 * @see dirichlet/support.h
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "support_collector.h"
#include "dirichlet/support.h"

namespace fem::bc {

/**
 * Constructs an empty named Dirichlet collector.
 *
 * @param name Collector name used by parser and loadcase selection.
 */
SupportCollector::SupportCollector(const std::string& name)
    : model::Collection<Dirichlet::Ptr>(name) {}

/**
 * Generates structural support equations from this collector.
 *
 * This compatibility overload intentionally delegates to the typed retrieval
 * path and selects only `Support` entries. Other Dirichlet definitions, such as
 * prescribed temperature, remain in the collector but do not enter the
 * structural constraint system.
 *
 * @param model_data Compiled model data used by structural supports.
 * @return Structural support equations in collector insertion order.
 */
constraint::Equations SupportCollector::get_equations(model::ModelData& model_data) {
    // Preserve the previous structural API while sharing all filtering and
    // application behavior with the generic typed implementation.
    return get_equations<Support>(model_data);
}

} // namespace fem::bc
