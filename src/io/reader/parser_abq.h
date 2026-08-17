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
 * The Abaqus reader retains only syntax state that has no persistent FEMaster
 * model counterpart. Nodal `*TRANSFORM` definitions are remembered until loads
 * and supports are materialized, while active load and boundary-condition
 * records are propagated between Abaqus steps independently of the mechanical
 * solution state.
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

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace fem::io::reader {

/**
 * @brief Logical Abaqus concentrated-load definition retained between steps.
 *
 * The original target spelling is preserved deliberately: Abaqus identifies a
 * load by the same node/node-set specification and generalized degree of
 * freedom when a later `OP=MOD` record changes it. Expansion to individual
 * FEMaster nodes is deferred until the current step is executed.
 */
struct ParserAbqCLoad {
    std::string target;
    int         dof = 0;
    Precision   magnitude = Precision(0);
    std::string amplitude;
    int         modified_step = 0;
};

/**
 * @brief Logical Abaqus displacement boundary condition retained between steps.
 *
 * Ranged `*BOUNDARY` input is stored as one record per generalized degree of
 * freedom so later modifications replace exactly the affected target/DOF pair.
 */
struct ParserAbqBoundary {
    std::string target;
    int         dof = 0;
    Precision   magnitude = Precision(0);
    std::string amplitude;
    int         modified_step = 0;
};

/**
 * @brief Logical supported Abaqus element-based distributed-load definition.
 *
 * The current reader supports `GRAV`; the type string is nevertheless retained
 * as part of the Abaqus load identity used by `OP=MOD`.
 */
struct ParserAbqDLoad {
    std::string              target;
    std::string              type;
    Precision                magnitude = Precision(0);
    std::array<Precision, 3> direction{Precision(0), Precision(0), Precision(0)};
    std::string              amplitude;
    int                      modified_step = 0;
};

/**
 * @brief Logical supported Abaqus surface-load definition retained between steps.
 *
 * Pressure and general traction share the same history container but remain
 * distinct through `type`. Tractions additionally retain their reference
 * direction, optional orientation and follower setting until materialization.
 */
struct ParserAbqDSLoad {
    std::string              surface;
    std::string              type;
    Precision                magnitude = Precision(0);
    std::array<Precision, 3> direction{Precision(0), Precision(0), Precision(0)};
    std::string              amplitude;
    std::string              orientation;
    std::string              follower;
    int                      modified_step = 0;
};

/**
 * @brief Syntax state required while reading one Abaqus deck.
 *
 * FEMaster load cases remain mechanically independent: every supported Abaqus
 * step starts from FEMaster's reference configuration and does not inherit
 * displacement, stress, material, contact or other converged solver state.
 *
 * Load and boundary-condition definitions are intentionally different. Abaqus
 * `OP=MOD` propagates and updates the corresponding logical records, while
 * `OP=NEW` clears that complete record category before the step supplies its new
 * definitions. At `*END STEP` the active snapshot is converted to a fresh
 * FEMaster collector for that one independent load case.
 */
struct ParserAbqState {
    // Abaqus nodal transform assignment retained between topology and data passes
    std::unordered_map<ID, std::string> node_transforms;

    // Active Abaqus load/BC definitions propagated independently of solver state
    std::vector<ParserAbqCLoad>    cloads;
    std::vector<ParserAbqBoundary> boundaries;
    std::vector<ParserAbqDLoad>    dloads;
    std::vector<ParserAbqDSLoad>   dsloads;

    // OP consistency within the currently open step
    std::string cload_op;
    std::string boundary_op;
    std::string dload_op;
    std::string dsload_op;

    // Monotonic identifier used for generated step-local FEMaster collectors
    int next_step_index = 1;

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

    // FEMaster collectors materialized from the current active snapshot
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
 * Parts, assemblies and mechanical state propagation between steps remain
 * deliberately outside the supported subset.
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
