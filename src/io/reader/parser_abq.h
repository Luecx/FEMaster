/**
 * @file parser_abq.h
 * @brief Declares the Abaqus syntax specialization of the staged FEMaster reader.
 *
 * `ParserAbq` reuses the model-construction pipeline implemented by `Parser` and
 * replaces only the command registration and activation performed in each
 * parser stage. This keeps allocation, topology finalization, element-local
 * enumeration, shell-normal construction and result-writer setup identical for
 * native and Abaqus input decks.
 *
 * The supported Abaqus syntax is intentionally introduced incrementally. The
 * reader currently accepts only the small model-definition subset registered in
 * `parser_abq.cpp`; unsupported Abaqus keywords remain hard parser errors.
 *
 * @see Parser
 * @see commands_abq::register_element
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include "parser.h"

namespace fem::io::reader {

/**
 * @brief Abaqus syntax specialization of the common staged input parser.
 *
 * The class does not own a second model-construction workflow. `Parser::run()`
 * remains responsible for the four dependency-ordered passes and all model-side
 * finalization between them. This specialization only defines which commands
 * are known in each pass and which of those commands are active.
 *
 * At the current implementation stage `HEADING`, `NODE`, `NSET` and `ELSET`
 * reuse the native FEMaster command definitions. `ELEMENT` is registered through
 * the Abaqus-specific command implementation because Abaqus element labels need
 * an explicit mapping onto FEMaster element formulations.
 *
 * The class contains no Abaqus part, assembly, step, load or material state yet;
 * those concepts are deliberately outside the current supported syntax.
 */
class ParserAbq : public Parser {
public:
    // Construction and destruction
    ParserAbq() = default;
    ~ParserAbq() override = default;

protected:
    // Stage-specific Abaqus command registration. The common Parser executes
    // these stages in the same order used for native FEMaster decks.
    void configure_count_stage   (io::dsl::Registry& registry, CountData& count) override;
    void configure_topology_stage(io::dsl::Registry& registry) override;
    void configure_field_stage   (io::dsl::Registry& registry) override;
    void configure_data_stage    (io::dsl::Registry& registry) override;
};

} // namespace fem::io::reader
