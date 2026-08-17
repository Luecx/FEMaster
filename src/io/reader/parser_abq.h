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
 * The Abaqus reader also owns the small amount of syntax state that has no
 * persistent FEMaster model counterpart. In particular, `*TRANSFORM` associates
 * coordinate systems with nodes for later load/support creation and `*STEP`
 * holds the currently parsed analysis procedure and its generated collector
 * names until `*END STEP` executes the corresponding FEMaster load case.
 *
 * @see Parser
 * @see commands::register_elastic
 * @see commands_abq::register_element
 * @see commands_abq::register_surface
 * @see commands_abq::register_orientation
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include "parser.h"

#include <string>
#include <unordered_map>

namespace fem::io::reader {

/**
 * @brief Transient syntax state required while reading one Abaqus deck.
 *
 * FEMaster stores coordinate systems on individual load/support objects rather
 * than attaching a transformed degree-of-freedom basis to nodes. The
 * `node_transforms` map therefore retains the Abaqus `*TRANSFORM` association
 * until step-local `*CLOAD` and `*BOUNDARY` records are converted to FEMaster
 * objects.
 *
 * The remaining members describe the currently open Abaqus step. Every step is
 * translated to one independent FEMaster load case. Loads and supports created
 * inside the step are collected under generated names and attached to that load
 * case when `*END STEP` is reached. Abaqus' persistent NLGEOM switch is retained
 * separately from the current-step flag because enabling geometric nonlinearity
 * affects all subsequent steps in an Abaqus/Standard analysis.
 */
struct ParserAbqState {
    // Abaqus nodal transform assignment retained between topology and data passes
    std::unordered_map<ID, std::string> node_transforms;

    // Deck-persistent syntax state
    int  next_step_index = 1;
    bool nlgeom_active   = false;

    // Currently open Abaqus step
    bool        step_active = false;
    int         step_index  = 0;
    int         max_increments = 100;
    bool        nlgeom       = false;
    bool        perturbation = false;
    Precision   step_period = Precision(1);
    std::string step_name;
    std::string step_amplitude;
    std::string procedure;

    // FEMaster collectors generated for the current step
    std::string load_collector;
    std::string support_collector;
    bool        load_collector_used    = false;
    bool        support_collector_used = false;
};

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
 * `DENSITY`, `ELASTIC` and Neo-Hooke `HYPERELASTIC`. Abaqus-specific
 * registrations handle element labels, surfaces, coordinate/orientation input,
 * thermal expansion, homogeneous sections and supported analysis/history data.
 *
 * Parts, assemblies and unsupported Abaqus procedures remain deliberately
 * outside the current subset.
 */
class ParserAbq : public Parser {
private:
    // Syntax-only state shared by Abaqus command callbacks in different passes
    ParserAbqState m_abq_state;

public:
    // Construction and destruction
    ParserAbq() = default;
    ~ParserAbq() override = default;

    // Mutable syntax state used by the dedicated Abaqus command registrations
    ParserAbqState& abaqus_state();
    const ParserAbqState& abaqus_state() const;

protected:
    // Stage-specific Abaqus command registration. The common Parser executes
    // these stages in the same order used for native FEMaster decks.
    void configure_count_stage   (io::dsl::Registry& registry, CountData& count) override;
    void configure_topology_stage(io::dsl::Registry& registry) override;
    void configure_field_stage   (io::dsl::Registry& registry) override;
    void configure_data_stage    (io::dsl::Registry& registry) override;
};

} // namespace fem::io::reader
