/**
 * @file parser_abq.cpp
 * @brief Implements Abaqus-specific command registration for the common parser stages.
 *
 * The Abaqus reader reuses native FEMaster registrations whenever the accepted
 * syntax and resulting model operation are identical. `ELASTIC` now uses the
 * shared command because the native syntax follows the Abaqus type meanings.
 * Format-specific handlers remain only where the external representation differs,
 * currently for `ELEMENT`, `SURFACE`, `ORIENTATION`, `EXPANSION`, `SOLID SECTION`
 * and `SHELL SECTION`.
 *
 * The supported material subset consists of `MATERIAL`, constant `DENSITY`,
 * isotropic or orthotropic `ELASTIC`, Neo-Hooke `HYPERELASTIC`, and constant
 * isotropic `EXPANSION`. Basic coordinate-defined rectangular orientations and
 * homogeneous solid, truss and shell sections are constructed in the topology
 * pass before the common parser assigns sections and enumerates element-local
 * data.
 *
 * No separate parsing pipeline is implemented here. The four passes, model
 * allocation, section assignment, element-local enumeration, shell-normal
 * construction and final data stage remain owned by `Parser`.
 *
 * @see Parser
 * @see commands::register_elastic
 * @see commands_abq::register_element
 * @see commands_abq::register_surface
 * @see commands_abq::register_orientation
 * @see commands_abq::register_expansion
 * @see commands_abq::register_solid_section
 * @see commands_abq::register_shell_section
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
#include "commands_abq/register_element.inl"
#include "commands_abq/register_expansion.inl"
#include "commands_abq/register_orientation.inl"
#include "commands_abq/register_shell_section.inl"
#include "commands_abq/register_solid_section.inl"
#include "commands_abq/register_surface.inl"

#include <algorithm>

namespace fem::io::reader {

/**
 * Configures the allocation pass for the currently supported Abaqus syntax.
 *
 * `NODE` and `ELEMENT` are the only active commands because model capacities
 * depend on their highest identifiers. All other supported keywords are still
 * registered in consume mode so a complete supported deck can be read without
 * mutating the temporary placeholder model. Abaqus surfaces use generated
 * internal surface identifiers and therefore require no surface count pass.
 *
 * @param registry Stage-local command registry.
 * @param count Allocation counters updated from parsed node and element ids.
 */
void ParserAbq::configure_count_stage(io::dsl::Registry& registry, CountData& count) {
    // Register every keyword currently accepted by the Abaqus reader
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
    commands_abq::register_solid_section(registry, model());
    commands_abq::register_shell_section(registry, model());

    // Execute only commands required to determine model allocation capacities
    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("NODE"   , io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELEMENT", io::dsl::ActiveMode::Active);
}

/**
 * Configures the topology pass for the currently supported Abaqus syntax.
 *
 * Nodes, elements, sets and surfaces are constructed together with material,
 * orientation and section definitions required before the common
 * section-assignment and element-enumeration steps. Material child commands
 * execute while their active `MATERIAL` parent selects the target material.
 *
 * @param registry Stage-local command registry.
 */
void ParserAbq::configure_topology_stage(io::dsl::Registry& registry) {
    // Register the supported model-definition commands
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
    commands_abq::register_solid_section(registry, model());
    commands_abq::register_shell_section(registry, model());

    // Construct all model data required before section assignment and element
    // enumeration. HEADING remains a recognized consume-only no-op.
    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("NODE"        , io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELEMENT"     , io::dsl::ActiveMode::Active);
    registry.set_active_mode("NSET"        , io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELSET"       , io::dsl::ActiveMode::Active);
    registry.set_active_mode("SURFACE"     , io::dsl::ActiveMode::Active);
    registry.set_active_mode("MATERIAL"    , io::dsl::ActiveMode::Active);
    registry.set_active_mode("DENSITY"     , io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELASTIC"     , io::dsl::ActiveMode::Active);
    registry.set_active_mode("HYPERELASTIC", io::dsl::ActiveMode::Active);
    registry.set_active_mode("EXPANSION"   , io::dsl::ActiveMode::Active);
    registry.set_active_mode("ORIENTATION" , io::dsl::ActiveMode::Active);
    registry.set_active_mode("SOLIDSECTION", io::dsl::ActiveMode::Active);
    registry.set_active_mode("SHELLSECTION", io::dsl::ActiveMode::Active);
}

/**
 * Configures the field pass for the currently supported Abaqus syntax.
 *
 * No Abaqus field-like keyword is supported yet. All known commands are
 * registered only in consume mode so the common field pass retains its position
 * between element enumeration and shell-normal construction without reapplying
 * topology, material, orientation or section mutations.
 *
 * @param registry Stage-local command registry.
 */
void ParserAbq::configure_field_stage(io::dsl::Registry& registry) {
    // Register the supported syntax so the complete deck can be consumed
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
    commands_abq::register_solid_section(registry, model());
    commands_abq::register_shell_section(registry, model());

    // No currently supported Abaqus keyword executes in the field pass
    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
}

/**
 * Configures the final data pass for the currently supported Abaqus syntax.
 *
 * The currently supported Abaqus subset contains only model-definition data,
 * all of which has already executed in the topology pass. Every command is
 * therefore consumed without execution while the common data-stage writer setup
 * completes. Loads, boundary conditions and analysis steps remain unsupported.
 *
 * @param registry Stage-local command registry.
 */
void ParserAbq::configure_data_stage(io::dsl::Registry& registry) {
    // Register the supported syntax so the complete deck can be consumed
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
    commands_abq::register_solid_section(registry, model());
    commands_abq::register_shell_section(registry, model());

    // No currently supported Abaqus keyword executes in the final data pass
    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
}

} // namespace fem::io::reader
