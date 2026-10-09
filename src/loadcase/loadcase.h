/**
 * @file loadcase.h
 * @brief Declares the base class for load-case execution.
 *
 * Load cases contain analysis settings and execute the associated solution
 * procedure. Loads, supports and their active named groups belong to ModelData;
 * run() assembles exactly the model state present when the analysis executes.
 *
 * @see src/loadcase/loadcase.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../core/castable.h"
#include "../constraints/constraint_groups.h"
#include "../model/model.h"
#include "../io/writer/writers.h"
#include "../io/writer/output_request_handler.h"

#include <memory>
#include <string>

namespace fem {
namespace loadcase {

/**
 * @brief Common runtime interface for finite-element analyses.
 *
 * A load case stores the analysis settings collected from consecutive input
 * commands and executes one concrete solution procedure against a model. The
 * parser owns active load cases through @ref Ptr and supplies their sequential
 * identifier, result writer and model when `Parser::begin_loadcase()` opens the
 * corresponding definition scope.
 *
 * The identifier and dependency pointers are therefore not yet bound on a
 * freshly constructed concrete load case and become valid before any setting
 * command can observe it. Derived classes provide the canonical input-deck type
 * name and implement the actual analysis in @ref run.
 *
 * Conditions and collector-activated conditions belong exclusively to ModelData. Derived
 * analyses assemble that current state during run(), without storing collector
 * lists, direct-condition copies or input-history transitions.
 *
 * Castable provides mutable and const runtime access to concrete implementations.
 *
 * The writer and model pointers are non-owning. Their lifetime is controlled by
 * the parser and covers the complete execution of every active load case.
 */
struct LoadCase : public Castable {
    // Owning polymorphic load-case pointer
    using Ptr = std::unique_ptr<LoadCase>;

    // Parser-assigned identity and shared analysis dependencies
    ID                         id     = -1;
    io::writer::ResultWriters* writer = nullptr;
    model::Model*              model  = nullptr;

    // Per-step output requests, provided basis fields and lazy recovered fields
    io::writer::OutputRequestHandler output;

    // Diagnostic settings
    bool report_constraints = false;

    // Construction and destruction
    virtual ~LoadCase() = default;

    // Concrete analysis identity and execution
    virtual std::string type_name() const = 0;
    virtual void run() = 0;

    // Assigned load-case identifier
    ID get_id() const { return id; }

protected:
    // Log the equation origins and sizes collected from the current model state
    void report_constraint_groups(const constraint::ConstraintGroups& groups) const;
};
} // namespace loadcase
} // namespace fem
