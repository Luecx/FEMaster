/**
 * @file loadcase.cpp
 * @brief Implements the base load-case utilities.
 *
 * Provides diagnostic helpers that report the assembled constraints for a
 * particular load case.
 *
 * @see src/loadcase/loadcase.h
 * @see src/model/constraint_groups.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "loadcase.h"

#include "../constraints/constraint_groups.h"

namespace fem {
namespace loadcase {

/**
 * @copydoc LoadCase::LoadCase
 */
LoadCase::LoadCase(ID case_id, reader::Writer* writer_in, model::Model* model_in)
    : id(case_id)
    , writer(writer_in)
    , model(model_in) {}

/**
 * @copydoc LoadCase::report_constraint_groups
 */
void LoadCase::report_constraint_groups(const constraint::ConstraintGroups& groups) const {
    if (!report_constraints) {
        return;
    }
    groups.report();
}

} // namespace loadcase
} // namespace fem
