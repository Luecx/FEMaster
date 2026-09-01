/**
 * @file model_robin.cpp
 * @brief Implements model-side registration of Robin boundary conditions.
 *
 * Robin conditions are stored in the same active load collector infrastructure
 * as Neumann conditions so input steps can select one named load set containing
 * heterogeneous natural and mixed boundary conditions. Registration happens
 * only after model compilation because boundary-condition regions address the
 * compiled assembly topology.
 *
 * The actual distinction between Neumann and Robin mathematical behavior is made
 * later by thermal load assembly, where Robin entries receive both the scalar
 * right-hand side and operator-triplet context.
 *
 * @see Model::add_load
 * @see bc::Robin
 * @see model_thermal_load.cpp
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "model.h"

#include "../core/logging.h"

namespace fem::model {

/**
 * Registers a Robin boundary condition in the currently active load collector.
 *
 * The model must already be compiled so region references are in the global
 * assembly identifier space. The function preserves polymorphic shared ownership
 * and does not assemble either side of the Robin condition; numerical evaluation
 * remains deferred until a thermal analysis requests the selected load collector.
 *
 * @param load Mixed boundary condition to append to the active collector.
 */
void Model::add_load(bc::Robin::Ptr load) {
    // Validate the compiled assembly state because Robin targets refer to
    // compiled node/surface topology rather than part-local semantic identifiers.
    logging::error(_data != nullptr && _data->compiled,
        "Model: loads require a compiled model");

    // Reject null load definitions before transferring shared ownership into the
    // collector.
    logging::error(load != nullptr,
        "Model: cannot add a null Robin load");

    // Boundary conditions are added to the collector currently selected by the
    // parser or direct model-building code. A missing active collector is an
    // invalid model-construction state rather than a silent no-op.
    logging::error(_data->load_cols.has_any() && _data->load_cols.get() != nullptr,
        "Model: no load collector is active");

    // Preserve the concrete Robin object through the common load collector. The
    // thermal assembly path later dispatches it through `bc::Robin`.
    _data->load_cols.get()->add(std::move(load));
}

} // namespace fem::model
