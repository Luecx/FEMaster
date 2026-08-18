/**
 * @file parser_abq.h
 * @brief Declares the staged Abaqus input reader and its parser-local state.
 *
 * `ParserAbq` specializes the command registration used by the common
 * `Parser::run()` pipeline for Abaqus input syntax. The class stores parser-only
 * information that is required while translating Abaqus model and history data
 * to FEMaster objects, including nodal `*TRANSFORM` assignments, active load and
 * boundary-condition definitions, and the currently open analysis step.
 *
 * @see Parser
 * @see commands_abq::register_transform
 * @see commands_abq::register_step
 * @see commands_abq::register_history
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
 * @brief Abaqus concentrated-load record used for step-wise load history.
 *
 * The target retains the original node or node-set specification and `dof`
 * identifies the generalized load component. `magnitude` is the active total
 * value. `previous_magnitude` stores the preceding total while a redefinition is
 * being materialized with the default `RAMP` behavior. `amplitude` references an
 * optional named amplitude and `modified_step` identifies the step in which the
 * active definition was last supplied.
 */
struct ParserAbqCLoad {
    std::string target;
    int         dof = 0;
    Precision   magnitude = Precision(0);
    Precision   previous_magnitude = Precision(0);
    std::string amplitude;
    int         modified_step = 0;
};

/**
 * @brief Abaqus displacement boundary-condition record used across steps.
 *
 * Each record represents one target and one generalized degree of freedom.
 * Ranged `*BOUNDARY` input is split into individual records so later `OP=MOD`
 * definitions can address each constrained DOF independently.
 */
struct ParserAbqBoundary {
    std::string target;
    int         dof = 0;
    Precision   magnitude = Precision(0);
    std::string amplitude;
    int         modified_step = 0;
};

/**
 * @brief Abaqus element-based distributed-load record used across steps.
 *
 * The supported record type is `GRAV`. Target, load type and direction form the
 * logical identity used for `OP=MOD`. `previous_magnitude` stores the preceding
 * total while a default transient `RAMP` is materialized.
 */
struct ParserAbqDLoad {
    std::string              target;
    std::string              type;
    Precision                magnitude = Precision(0);
    Precision                previous_magnitude = Precision(0);
    std::array<Precision, 3> direction{Precision(0), Precision(0), Precision(0)};
    std::string              amplitude;
    int                      modified_step = 0;
};

/**
 * @brief Abaqus surface-load record used across steps.
 *
 * Pressure and general traction records store their surface, load type and
 * scalar magnitude. General traction additionally stores its direction,
 * orientation and follower setting. `previous_magnitude` stores the preceding
 * total while a default transient `RAMP` is materialized.
 */
struct ParserAbqDSLoad {
    std::string              surface;
    std::string              type;
    Precision                magnitude = Precision(0);
    Precision                previous_magnitude = Precision(0);
    std::array<Precision, 3> direction{Precision(0), Precision(0), Precision(0)};
    std::string              amplitude;
    std::string              orientation;
    std::string              follower;
    int                      modified_step = 0;
};

/**
 * @brief Parser-local state for one Abaqus input deck.
 *
 * The state stores nodal transformation assignments, active load and boundary
 * definitions, `OP` selections for the current step, and the parameters required
 * to construct the active FEMaster load case. Load and boundary definitions are
 * propagated between steps, whereas each FEMaster load case uses its own
 * mechanical initial state.
 */
struct ParserAbqState {
    // Nodal coordinate systems assigned by Abaqus *TRANSFORM
    std::unordered_map<ID, std::string> node_transforms;

    // Active Abaqus load and boundary-condition definitions
    std::vector<ParserAbqCLoad>    cloads;
    std::vector<ParserAbqBoundary> boundaries;
    std::vector<ParserAbqDLoad>    dloads;
    std::vector<ParserAbqDSLoad>   dsloads;

    // OP values used by each definition category in the open step
    std::string cload_op;
    std::string boundary_op;
    std::string dload_op;
    std::string dsload_op;

    // Monotonic identifier for generated step-local FEMaster collectors
    int next_step_index = 1;

    // Active Abaqus step and procedure parameters
    bool        step_active = false;
    int         step_index  = 0;
    int         max_increments = 100;
    bool        nlgeom       = false;
    bool        perturbation = false;
    Precision   step_period = Precision(1);
    std::string step_name;
    std::string step_amplitude;
    std::string procedure;

    // FEMaster collector names created for the active step snapshot
    std::string load_collector;
    std::string support_collector;
    bool        load_collector_used    = false;
    bool        support_collector_used = false;
};

/**
 * @brief Abaqus specialization of the common staged FEMaster input parser.
 *
 * The class registers Abaqus-compatible commands for allocation, topology,
 * field, and analysis-data stages while delegating the stage execution and model
 * finalization to `Parser`. Native command registrations are reused where their
 * syntax and semantics match the supported Abaqus subset.
 */
class ParserAbq : public Parser {
private:
    // Parser-local Abaqus syntax and history state
    ParserAbqState m_abq_state;

public:
    // Construction and destruction
    ParserAbq() = default;
    ~ParserAbq() override = default;

    // Abaqus state shared by format-specific command registrations
    ParserAbqState& abaqus_state();
    const ParserAbqState& abaqus_state() const;

protected:
    // Stage-specific command registration
    void configure_count_stage   (io::dsl::Registry& registry, CountData& count) override;
    void configure_topology_stage(io::dsl::Registry& registry) override;
    void configure_field_stage   (io::dsl::Registry& registry) override;
    void configure_data_stage    (io::dsl::Registry& registry) override;
};

} // namespace fem::io::reader
