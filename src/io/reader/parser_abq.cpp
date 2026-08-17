/**
 * @file parser_abq.cpp
 * @brief Implements Abaqus-specific stage command registration.
 *
 * The initial reader supports only HEADING, NODE, ELEMENT, NSET and ELSET.
 * HEADING, NODE, NSET and ELSET reuse the native FEMaster command definitions;
 * ELEMENT has a dedicated Abaqus registration for Abaqus element type names.
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#include "parser_abq.h"

#include "commands/register_node_count.inl"
#include "commands/register_heading.inl"
#include "commands/register_node.inl"
#include "commands/register_nset.inl"
#include "commands/register_elset.inl"
#include "commands_abq/register_element.inl"

#include <algorithm>

namespace fem::io::reader {

void ParserAbq::configure_count_stage(io::dsl::Registry& registry, CountData& count) {
    commands::register_heading(registry);
    commands::register_node_count(registry, [&count](ID id) {
        count.highest_node_id = std::max(count.highest_node_id, static_cast<int>(id));
    });
    commands_abq::register_element_count(registry, [&count](ID id) {
        count.highest_element_id = std::max(count.highest_element_id, static_cast<int>(id));
    });
    commands::register_nset(registry, model());
    commands::register_elset(registry, model());

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("NODE", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELEMENT", io::dsl::ActiveMode::Active);
}

void ParserAbq::configure_topology_stage(io::dsl::Registry& registry) {
    commands::register_heading(registry);
    commands::register_node(registry, model());
    commands_abq::register_element(registry, model());
    commands::register_nset(registry, model());
    commands::register_elset(registry, model());

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("NODE", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELEMENT", io::dsl::ActiveMode::Active);
    registry.set_active_mode("NSET", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELSET", io::dsl::ActiveMode::Active);
}

void ParserAbq::configure_field_stage(io::dsl::Registry& registry) {
    commands::register_heading(registry);
    commands::register_node(registry, model());
    commands_abq::register_element(registry, model());
    commands::register_nset(registry, model());
    commands::register_elset(registry, model());

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
}

void ParserAbq::configure_data_stage(io::dsl::Registry& registry) {
    commands::register_heading(registry);
    commands::register_node(registry, model());
    commands_abq::register_element(registry, model());
    commands::register_nset(registry, model());
    commands::register_elset(registry, model());

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
}

} // namespace fem::io::reader
