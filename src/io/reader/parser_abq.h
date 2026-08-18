/**
 * @file parser_abq.h
 * @brief Declares the staged Abaqus input reader and its parser-local state.
 *
 * `ParserAbq` specializes the command registration used by the common
 * `Parser::run()` pipeline for the supported Abaqus input syntax. The reader
 * accepts at most one analysis `*STEP`, so no mechanical or load-definition
 * history is stored between steps.
 *
 * Parser-local state is limited to nodal `*TRANSFORM` assignments and parameters
 * required while translating the active analysis step to a FEMaster load case.
 *
 * @see Parser
 * @see commands_abq::register_transform
 * @see commands_abq::register_step
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include "parser.h"

#include <string>
#include <unordered_map>
#include <utility>

namespace fem::io::reader {

/**
 * @brief Parser-local state for one supported Abaqus input deck.
 *
 * The state stores nodal coordinate transformations and the control parameters
 * of the single supported analysis step. Loads and boundary conditions are
 * translated directly into step-local FEMaster collectors while their input
 * cards are parsed.
 */
struct ParserAbqState {
    // Nodal coordinate systems assigned by Abaqus *TRANSFORM
    std::unordered_map<ID, std::string> node_transforms;

    // Single supported Abaqus analysis step
    bool        step_seen      = false;
    bool        step_active    = false;
    int         max_increments = 100;
    bool        nlgeom         = false;
    bool        perturbation   = false;
    Precision   step_period    = Precision(1);
    std::string step_amplitude;
    std::string procedure;
};

/**
 * @brief Abaqus specialization of the common staged FEMaster input parser.
 *
 * The class registers Abaqus-compatible commands for allocation, topology,
 * field, and analysis-data stages while delegating stage execution and model
 * finalization to `Parser`. The supported analysis subset contains at most one
 * `*STEP`; step-local loads and supports are written directly to FEMaster
 * collectors.
 */
class ParserAbq : public Parser {
private:
    // Parser-local Abaqus syntax state
    ParserAbqState m_abq_state;

public:
    // Construction and destruction
    ParserAbq() = default;
    ~ParserAbq() override = default;

    // Abaqus state shared by format-specific command registrations
    ParserAbqState& abaqus_state();
    const ParserAbqState& abaqus_state() const;

    // Map a load amplitude to the representation used by the active load case
    std::pair<Precision, std::string> resolve_load_amplitude(const std::string& amplitude);

protected:
    // Stage-specific command registration
    void configure_count_stage   (io::dsl::Registry& registry, CountData& count) override;
    void configure_topology_stage(io::dsl::Registry& registry) override;
    void configure_field_stage   (io::dsl::Registry& registry) override;
    void configure_data_stage    (io::dsl::Registry& registry) override;
};

} // namespace fem::io::reader
