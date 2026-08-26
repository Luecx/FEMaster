/**
 * @file parser_abq.h
 * @brief Declares the Abaqus reader and its parser-local analysis state.
 *
 * `ParserAbq` registers the supported Abaqus grammar once and then processes the
 * resulting `dsl::Deck` in an explicit dependency order. Materials and Part
 * topology are created first, `Model::compile()` forms the dense assembly, and
 * assembly regions/transforms plus the final analysis step execute afterwards.
 *
 * Parser-local state is limited to nodal `*TRANSFORM` assignments and controls
 * required while translating one supported analysis step.
 *
 * @see Parser
 * @see commands_abq::register_transform
 * @see commands_abq::register_step
 *
 * @author Finn Eggers
 * @date 26.08.2026
 */

#pragma once

#include "parser.h"

#include <string>
#include <unordered_map>
#include <utility>

namespace fem::io::reader {

/**
 * @brief Runtime state required while translating one supported Abaqus step.
 */
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

/**
 * @brief Abaqus dialect specialization of the common parsed-deck reader.
 */
class ParserAbq : public Parser {
private:
    ParserAbqState m_abq_state;

public:
    ParserAbq();
    ~ParserAbq() override = default;

    ParserAbqState& abaqus_state();
    const ParserAbqState& abaqus_state() const;

    std::pair<Precision, std::string> resolve_load_amplitude(const std::string& amplitude);

protected:
    void register_commands(io::dsl::Registry& registry) override;
    void process_deck(const io::dsl::Deck&                  deck,
                      const std::string&                    input_path,
                      const std::string&                    output_path,
                      const io::writer::WriterFileFormats& writer_formats) override;

private:
    void register_common_commands(io::dsl::Registry& registry);
};

} // namespace fem::io::reader
