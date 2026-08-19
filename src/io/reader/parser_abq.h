/**
 * @file parser_abq.h
 * @brief Declares the Abaqus input reader and its parser-local state.
 *
 * `ParserAbq` specializes the semantic parser passes used by the common
 * `Parser::run()` pipeline. Global material resources are collected before
 * topology, nodes and elements are constructed in reusable parts, and the
 * common `Model::compile()` boundary creates dense solver data before assembly
 * sets/surfaces, transforms and analysis commands execute.
 *
 * The reader accepts at most one analysis `*STEP`, so no mechanical or
 * load-definition history is stored between steps. Parser-local state is limited
 * to nodal `*TRANSFORM` assignments and controls required while translating that
 * step to a FEMaster load case.
 *
 * @see Parser
 * @see commands_abq::register_transform
 * @see commands_abq::register_step
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "parser.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

namespace fem::io::reader {

struct ParserAbqState {
    std::unordered_map<ID, std::string> node_transforms;

    bool        step_seen      = false;
    bool        step_active    = false;
    int         max_increments = 100;
    bool        nlgeom         = false;
    bool        perturbation   = false;
    Precision   step_period    = Precision(1);
    std::string step_amplitude;
};

class ParserAbq : public Parser {
private:
    ParserAbqState m_abq_state;

public:
    ParserAbq() = default;
    ~ParserAbq() override = default;

    ParserAbqState& abaqus_state();
    const ParserAbqState& abaqus_state() const;

    std::pair<Precision, std::string> resolve_load_amplitude(const std::string& amplitude);

protected:
    void configure_definition_pass(io::dsl::Registry& registry) override;
    void configure_topology_pass  (io::dsl::Registry& registry) override;
    void configure_assembly_pass  (io::dsl::Registry& registry) override;
    void configure_analysis_pass  (io::dsl::Registry& registry) override;

private:
    void register_common_commands(io::dsl::Registry& registry,
                                  std::shared_ptr<bool> assembly_scope);
};

} // namespace fem::io::reader
