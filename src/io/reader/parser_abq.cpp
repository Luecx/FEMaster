/**
 * @file parser_abq.cpp
 * @brief Implements Abaqus-specific command registration for the common parser stages.
 *
 * The Abaqus reader reuses native FEMaster registrations whenever the accepted
 * syntax and resulting model operation are identical. Format-specific handlers
 * remain only where the external representation or Abaqus history semantics
 * differ.
 *
 * Model-definition data are constructed in the topology pass. This includes
 * nodes/elements/sets/surfaces, materials and sections, tabular amplitudes,
 * material orientations and nodal `*TRANSFORM` coordinate systems. Analysis
 * steps and their load/support histories execute only in the final data pass,
 * after the common parser has assigned sections and initialized all element-local
 * data.
 *
 * Supported steps remain mechanically independent FEMaster load cases. Abaqus
 * load and boundary-condition definitions nevertheless propagate logically
 * between steps and are materialized as a fresh collector snapshot at each
 * `*END STEP`, so `OP=MOD` and `OP=NEW` do not require solver-state propagation.
 *
 * @see Parser
 * @see ParserAbqState
 * @see commands_abq::register_step
 * @see commands_abq::register_history
 * @see commands_abq::register_transform
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#include "parser_abq.h"

#include "commands/register_density.inl"
#include "commands/register_elastic.inl"
#include "commands/register_elset.inl"
#include "commands/register_heading.inl"
#include "commands/register_hyperelastic.inl"
#include "commands/register_material.inl"
#include "commands/register_node.inl"
#include "commands/register_node_count.inl"
#include "commands/register_nset.inl"
#include "commands_abq/register_amplitude.inl"
#include "commands_abq/register_boundary.inl"
#include "commands_abq/register_cload.inl"
#include "commands_abq/register_dload.inl"
#include "commands_abq/register_dsload.inl"
#include "commands_abq/register_element.inl"
#include "commands_abq/register_expansion.inl"
#include "commands_abq/register_history.inl"
#include "commands_abq/register_orientation.inl"
#include "commands_abq/register_shell_section.inl"
#include "commands_abq/register_solid_section.inl"
#include "commands_abq/register_step.inl"
#include "commands_abq/register_surface.inl"
#include "commands_abq/register_transform.inl"

#include <algorithm>

namespace fem::io::reader {

ParserAbqState& ParserAbq::abaqus_state() {
    return m_abq_state;
}

const ParserAbqState& ParserAbq::abaqus_state() const {
    return m_abq_state;
}

/**
 * Configures the allocation pass for the supported Abaqus syntax.
 *
 * `NODE` and `ELEMENT` are the only active commands because model capacities
 * depend on their highest identifiers. Every other supported keyword is still
 * registered in consume mode so complete supported decks can be traversed
 * without mutating the placeholder count-stage model.
 *
 * @param registry Stage-local command registry.
 * @param count Allocation counters updated from parsed node and element ids.
 */
void ParserAbq::configure_count_stage(io::dsl::Registry& registry, CountData& count) {
    // Model-definition syntax
    commands::register_heading(registry);
    commands::register_node_count(registry, [&count](ID id) {
        count.highest_node_id = std::max(count.highest_node_id, static_cast<int>(id));
    });
    commands_abq::register_element_count(registry, [&count](ID id) {
        count.highest_element_id = std::max(count.highest_element_id, static_cast<int>(id));
    });
    commands::register_nset(registry, model());
    commands::register_elset(registry, model());
    commands_abq::register_surface(registry, model());
    commands::register_material(registry, model());
    commands::register_density(registry, model());
    commands::register_elastic(registry, model());
    commands::register_hyperelastic(registry, model());
    commands_abq::register_expansion(registry, model());
    commands_abq::register_orientation(registry, model());
    commands_abq::register_transform(registry, *this);
    commands_abq::register_amplitude(registry, model());
    commands_abq::register_solid_section(registry, model());
    commands_abq::register_shell_section(registry, model());

    // Analysis/history syntax is recognized but inert during allocation
    commands_abq::register_step(registry, *this);
    commands_abq::register_cload(registry, *this);
    commands_abq::register_boundary(registry, *this);
    commands_abq::register_dload(registry, *this);
    commands_abq::register_dsload(registry, *this);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("NODE"   , io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELEMENT", io::dsl::ActiveMode::Active);
}

/**
 * Configures the topology pass for the supported Abaqus syntax.
 *
 * All persistent model data required by later load cases are built here,
 * including amplitudes and nodal transforms. Abaqus step/history commands remain
 * consume-only until the final pass so no solver or logical history state is
 * changed before model topology, sections and element-local enumeration are
 * complete.
 *
 * @param registry Stage-local command registry.
 */
void ParserAbq::configure_topology_stage(io::dsl::Registry& registry) {
    // A ParserAbq instance may parse more than one deck during its lifetime.
    // Reset all format-only state before constructing the new model.
    m_abq_state = ParserAbqState{};

    // Persistent model-definition syntax
    commands::register_heading(registry);
    commands::register_node(registry, model());
    commands_abq::register_element(registry, model());
    commands::register_nset(registry, model());
    commands::register_elset(registry, model());
    commands_abq::register_surface(registry, model());
    commands::register_material(registry, model());
    commands::register_density(registry, model());
    commands::register_elastic(registry, model());
    commands::register_hyperelastic(registry, model());
    commands_abq::register_expansion(registry, model());
    commands_abq::register_orientation(registry, model());
    commands_abq::register_transform(registry, *this);
    commands_abq::register_amplitude(registry, model());
    commands_abq::register_solid_section(registry, model());
    commands_abq::register_shell_section(registry, model());

    // History syntax must still be registered so the topology pass can traverse
    // complete Abaqus decks without executing a load case.
    commands_abq::register_step(registry, *this);
    commands_abq::register_cload(registry, *this);
    commands_abq::register_boundary(registry, *this);
    commands_abq::register_dload(registry, *this);
    commands_abq::register_dsload(registry, *this);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("NODE"          , io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELEMENT"       , io::dsl::ActiveMode::Active);
    registry.set_active_mode("NSET"          , io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELSET"         , io::dsl::ActiveMode::Active);
    registry.set_active_mode("SURFACE"       , io::dsl::ActiveMode::Active);
    registry.set_active_mode("MATERIAL"      , io::dsl::ActiveMode::Active);
    registry.set_active_mode("DENSITY"       , io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELASTIC"       , io::dsl::ActiveMode::Active);
    registry.set_active_mode("HYPERELASTIC"  , io::dsl::ActiveMode::Active);
    registry.set_active_mode("EXPANSION"     , io::dsl::ActiveMode::Active);
    registry.set_active_mode("ORIENTATION"   , io::dsl::ActiveMode::Active);
    registry.set_active_mode("TRANSFORM"     , io::dsl::ActiveMode::Active);
    registry.set_active_mode("AMPLITUDE"     , io::dsl::ActiveMode::Active);
    registry.set_active_mode("SOLIDSECTION"  , io::dsl::ActiveMode::Active);
    registry.set_active_mode("SHELLSECTION"  , io::dsl::ActiveMode::Active);
}

/**
 * Configures the field pass for the supported Abaqus syntax.
 *
 * No currently supported Abaqus keyword creates an enumeration-dependent field.
 * Every known command is registered only in consume mode so the common field
 * pass retains its dependency position without replaying topology or history
 * mutations.
 *
 * @param registry Stage-local command registry.
 */
void ParserAbq::configure_field_stage(io::dsl::Registry& registry) {
    commands::register_heading(registry);
    commands::register_node(registry, model());
    commands_abq::register_element(registry, model());
    commands::register_nset(registry, model());
    commands::register_elset(registry, model());
    commands_abq::register_surface(registry, model());
    commands::register_material(registry, model());
    commands::register_density(registry, model());
    commands::register_elastic(registry, model());
    commands::register_hyperelastic(registry, model());
    commands_abq::register_expansion(registry, model());
    commands_abq::register_orientation(registry, model());
    commands_abq::register_transform(registry, *this);
    commands_abq::register_amplitude(registry, model());
    commands_abq::register_solid_section(registry, model());
    commands_abq::register_shell_section(registry, model());
    commands_abq::register_step(registry, *this);
    commands_abq::register_cload(registry, *this);
    commands_abq::register_boundary(registry, *this);
    commands_abq::register_dload(registry, *this);
    commands_abq::register_dsload(registry, *this);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
}

/**
 * Configures the final Abaqus analysis/history pass.
 *
 * Persistent model-definition keywords are consumed without execution because
 * they have already been applied in topology. Step/procedure cards update the
 * logical Abaqus load/BC history; `register_history` then materializes the active
 * snapshot into one fresh FEMaster collector at `*END STEP` and executes the
 * mechanically independent load case.
 *
 * @param registry Stage-local command registry.
 */
void ParserAbq::configure_data_stage(io::dsl::Registry& registry) {
    // Persistent model syntax remains recognized but inert
    commands::register_heading(registry);
    commands::register_node(registry, model());
    commands_abq::register_element(registry, model());
    commands::register_nset(registry, model());
    commands::register_elset(registry, model());
    commands_abq::register_surface(registry, model());
    commands::register_material(registry, model());
    commands::register_density(registry, model());
    commands::register_elastic(registry, model());
    commands::register_hyperelastic(registry, model());
    commands_abq::register_expansion(registry, model());
    commands_abq::register_orientation(registry, model());
    commands_abq::register_transform(registry, *this);
    commands_abq::register_amplitude(registry, model());
    commands_abq::register_solid_section(registry, model());
    commands_abq::register_shell_section(registry, model());

    // Analysis/history syntax executes only after model finalization
    commands_abq::register_step(registry, *this);
    commands_abq::register_cload(registry, *this);
    commands_abq::register_boundary(registry, *this);
    commands_abq::register_dload(registry, *this);
    commands_abq::register_dsload(registry, *this);

    // Replace only STEP entry/ENDSTEP finalization semantics. Procedure variants
    // from register_step remain unchanged.
    commands_abq::register_history(registry, *this);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("STEP"               , io::dsl::ActiveMode::Active);
    registry.set_active_mode("STATIC"             , io::dsl::ActiveMode::Active);
    registry.set_active_mode("FREQUENCY"          , io::dsl::ActiveMode::Active);
    registry.set_active_mode("BUCKLE"             , io::dsl::ActiveMode::Active);
    registry.set_active_mode("DYNAMIC"            , io::dsl::ActiveMode::Active);
    registry.set_active_mode("STEADYSTATEDYNAMICS", io::dsl::ActiveMode::Active);
    registry.set_active_mode("CLOAD"              , io::dsl::ActiveMode::Active);
    registry.set_active_mode("BOUNDARY"           , io::dsl::ActiveMode::Active);
    registry.set_active_mode("DLOAD"              , io::dsl::ActiveMode::Active);
    registry.set_active_mode("DSLOAD"             , io::dsl::ActiveMode::Active);
    registry.set_active_mode("ENDSTEP"            , io::dsl::ActiveMode::Active);
}

} // namespace fem::io::reader
