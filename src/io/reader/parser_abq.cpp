/**
 * @file parser_abq.cpp
 * @brief Implements Abaqus-specific command registration for the common parser stages.
 *
 * The Abaqus reader currently recognizes only `HEADING`, `NODE`, `ELEMENT`,
 * `NSET` and `ELSET`. `HEADING`, `NODE`, `NSET` and `ELSET` reuse the native
 * FEMaster command definitions without format-specific changes. `ELEMENT` uses
 * the dedicated Abaqus registration because external element labels are mapped
 * onto FEMaster element formulations.
 *
 * No separate parsing pipeline is implemented here. The four passes, model
 * allocation, section assignment, element-local enumeration, shell-normal
 * construction and final data stage remain owned by `Parser`.
 *
 * @see Parser
 * @see commands_abq::register_element
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#include "parser_abq.h"

#include "commands/register_elset.inl"
#include "commands/register_heading.inl"
#include "commands/register_node.inl"
#include "commands/register_node_count.inl"
#include "commands/register_nset.inl"
#include "commands_abq/register_element.inl"

#include <algorithm>

namespace fem::io::reader {

/**
 * Configures the allocation pass for the currently supported Abaqus syntax.
 *
 * `NODE` and `ELEMENT` are the only active commands because model capacities
 * depend on their highest identifiers. `HEADING`, `NSET` and `ELSET` are still
 * registered so the DSL engine can consume a complete supported deck without
 * executing set mutations on the temporary placeholder model.
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

    // Execute only commands required to determine model allocation capacities
    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("NODE"   , io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELEMENT", io::dsl::ActiveMode::Active);
}

/**
 * Configures the topology pass for the currently supported Abaqus syntax.
 *
 * Nodes, elements and named node/element sets are constructed in the allocated
 * FEMaster model. `HEADING` remains a recognized no-op command. Element labels
 * are resolved through the Abaqus-specific registration before topology is
 * finalized by the common parser pipeline.
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

    // Construct topology and sets exactly once in the topology pass
    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("NODE"   , io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELEMENT", io::dsl::ActiveMode::Active);
    registry.set_active_mode("NSET"   , io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELSET"  , io::dsl::ActiveMode::Active);
}

/**
 * Configures the field pass for the currently supported Abaqus syntax.
 *
 * No Abaqus field-like keyword is supported yet. All currently known commands
 * are therefore registered only in consume mode so the common field pass keeps
 * its position between element enumeration and shell-normal construction
 * without reapplying topology mutations.
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

    // No currently supported Abaqus keyword executes in the field pass
    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
}

/**
 * Configures the final data pass for the currently supported Abaqus syntax.
 *
 * Loads, boundary conditions, materials, sections and analysis steps are not yet
 * part of the accepted Abaqus subset. The existing topology keywords are only
 * consumed so the common data-stage writer setup can complete without applying
 * model mutations a second time.
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

    // No currently supported Abaqus keyword executes in the final data pass
    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
}

} // namespace fem::io::reader
