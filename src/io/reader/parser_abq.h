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
 * reader currently covers basic topology and sets, element- and node-based
 * surfaces, core elastic material properties and homogeneous solid/truss/shell
 * sections. Unsupported Abaqus keywords remain hard parser errors.
 *
 * @see Parser
 * @see commands_abq::register_element
 * @see commands_abq::register_surface
 * @see commands_abq::register_elastic
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
 * Native registrations are reused where Abaqus input maps without semantic
 * changes, including `HEADING`, `NODE`, `NSET`, `ELSET`, `MATERIAL`, constant
 * `DENSITY` and Neo-Hooke `HYPERELASTIC`. Abaqus-specific registrations handle
 * element labels, surfaces, linear elasticity, thermal expansion and homogeneous
 * section syntax where the external representation differs from FEMaster.
 *
 * Parts, assemblies, analysis steps, loads and boundary conditions remain
 * deliberately outside the current supported subset.
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
