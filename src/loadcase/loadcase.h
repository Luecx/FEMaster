/******************************************************************************
 * @file loadcase.h
 * @brief Declares the base class for load-case execution.
 *
 * Load cases encapsulate the application of supports and loads on a model while
 * delegating solver specifics to derived classes.
 *
 * @see src/loadcase/loadcase.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#include "../model/model.h"
#include "../constraints/constraint_groups.h"
#include "../reader/writer.h"

namespace fem {
namespace loadcase {

/******************************************************************************
 * @struct LoadCase
 * @brief Abstract base representing a simulation load case.
 ******************************************************************************/
struct LoadCase {
    const ID id;              ///< Load-case identifier.
    reader::Writer* writer;   ///< Output writer associated with the load case.
    model::Model* model;      ///< Model on which the load case operates.
    bool report_constraints = false; ///< Flag to print constraint reports.

    /******************************************************************************
     * @brief Constructs a load case with identifiers and dependencies.
     *
     * @param case_id Identifier of the load case.
     * @param writer_in Writer used for logging or output.
     * @param model_in Model instance to operate on.
     ******************************************************************************/
    LoadCase(ID case_id, reader::Writer* writer_in, model::Model* model_in);

    /******************************************************************************
     * @brief Executes the load case.
     ******************************************************************************/
    virtual void run() = 0;

    /// Returns the load-case identifier.
    ID get_id() const { return id; }

protected:
    /******************************************************************************
     * @brief Logs a summary of the active constraint groups.
     *
     * @param groups Constraint groups extracted from the model.
     ******************************************************************************/
    void report_constraint_groups(const model::ConstraintGroups& groups) const;
};

} // namespace loadcase
} // namespace fem
