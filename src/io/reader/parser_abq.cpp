/**
 * @file parser_abq.cpp
 * @brief Implements staged command registration for Abaqus input decks.
 *
 * The Abaqus reader uses the common FEMaster parser lifecycle and configures the
 * commands that are recognized and executed in each stage. Topology data such as
 * nodes, elements, sets, materials, sections, amplitudes, orientations and nodal
 * transforms are created before analysis data. Step procedures and step-local
 * load or boundary definitions are evaluated in the final data stage.
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
 * Configures the allocation stage for the supported Abaqus syntax.
 *
 * Node and element identifiers are evaluated to determine model capacities. All
 * other supported commands are registered in consume-only mode so the complete
 * input structure can be traversed without creating model data in this stage.
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

    // Analysis syntax is recognized but not executed during allocation
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
 * Configures the topology stage for the supported Abaqus syntax.
 *
 * Persistent model data are created in this stage, including topology, sets,
 * materials, sections, amplitudes, material orientations and nodal transforms.
 * Analysis-step commands remain consume-only until the final data stage.
 *
 * @param registry Stage-local command registry.
 */
void ParserAbq::configure_topology_stage(io::dsl::Registry& registry) {
    // Start each deck with an empty parser-local Abaqus state.
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

    // Analysis syntax remains inactive while topology is constructed
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
 * Configures the field stage for the supported Abaqus syntax.
 *
 * The currently supported Abaqus subset does not define additional
 * enumeration-dependent fields. All known commands are therefore registered in
 * consume-only mode for this stage.
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
 * Configures the analysis-data stage for the supported Abaqus syntax.
 *
 * Model-definition commands are consumed without execution. Step procedures,
 * loads and boundary conditions are active; `register_history` materializes the
 * active logical definitions into FEMaster collectors when `*END STEP` is read.
 *
 * @param registry Stage-local command registry.
 */
void ParserAbq::configure_data_stage(io::dsl::Registry& registry) {
    // Persistent model syntax is recognized but inactive in this stage
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

    // Analysis and load-definition syntax
    commands_abq::register_step(registry, *this);
    commands_abq::register_cload(registry, *this);
    commands_abq::register_boundary(registry, *this);
    commands_abq::register_dload(registry, *this);
    commands_abq::register_dsload(registry, *this);
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
