/**
 * @file parser_abq.cpp
 * @brief Implements dependency-ordered Abaqus input-deck parsing.
 *
 * Abaqus decks are replayed in four semantic passes. Definitions are loaded
 * first, followed by part/instance topology. `Model::compile()` then flattens
 * that topology, including part-local point-mass features, exactly once.
 * Assembly sets, surfaces and nodal transforms are materialized against the
 * compiled assembly before the final step/load pass.
 *
 * Each pass registers the same command grammar but activates only commands whose
 * dependencies are satisfied at that point. All other keywords are consumed so
 * the original Abaqus nesting remains intact without executing callbacks too
 * early.
 *
 * @author Finn Eggers
 * @date 24.08.2026
 */

#include "parser_abq.h"

#include "commands/register_density.inl"
#include "commands/register_elastic.inl"
#include "commands/register_elset.inl"
#include "commands/register_end_part.inl"
#include "commands/register_heading.inl"
#include "commands/register_hyperelastic.inl"
#include "commands/register_material.inl"
#include "commands/register_node.inl"
#include "commands/register_nset.inl"
#include "commands/register_part.inl"
#include "commands/register_surface.inl"
#include "commands_abq/register_amplitude.inl"
#include "commands_abq/register_assembly.inl"
#include "commands_abq/register_boundary.inl"
#include "commands_abq/register_cload.inl"
#include "commands_abq/register_dload.inl"
#include "commands_abq/register_dsload.inl"
#include "commands_abq/register_element.inl"
#include "commands_abq/register_end_assembly.inl"
#include "commands_abq/register_end_instance.inl"
#include "commands_abq/register_equation.inl"
#include "commands_abq/register_expansion.inl"
#include "commands_abq/register_instance.inl"
#include "commands_abq/register_mass.inl"
#include "commands_abq/register_orientation.inl"
#include "commands_abq/register_shell_section.inl"
#include "commands_abq/register_solid_section.inl"
#include "commands_abq/register_step.inl"
#include "commands_abq/register_transform.inl"

#include "../dsl/registry.h"
#include "../../bc/amplitude.h"
#include "../../feature/point_mass.h"
#include "../../loadcase/loadcase.h"
#include "../../model/element/point.h"
#include "../../model/model.h"

#include <memory>
#include <string>
#include <utility>

namespace fem::io::reader {

ParserAbqState& ParserAbq::abaqus_state() { return m_abq_state; }
const ParserAbqState& ParserAbq::abaqus_state() const { return m_abq_state; }

/**
 * Creates one concentrated PointMass feature from a compiled Abaqus MASS element.
 *
 * This helper is used only for assembly-level `*MASS` records whose target
 * PointElement already lives in dense assembly space. The PointElement is used
 * only to resolve its single target node; the resulting PointMass remains a
 * purely nodal feature with no retained element association.
 *
 * @param element_id Dense compiled PointElement identifier.
 * @param mass Isotropic translational mass assigned by Abaqus `*MASS`.
 */
void ParserAbq::add_abaqus_mass_feature(ID element_id, Precision mass) {
    auto& data = *model()._data;

    logging::error(data.compiled,
        "MASS: point-mass features require a compiled model");
    logging::error(element_id >= 0 && static_cast<std::size_t>(element_id) < data.elements.size() && data.elements[static_cast<std::size_t>(element_id)] != nullptr,
        "MASS: compiled element ", element_id, " does not exist");

    auto* point = data.elements[static_cast<std::size_t>(element_id)]->as<model::PointElement>();
    logging::error(point != nullptr,
        "MASS: element ", element_id, " is not a MASS element");

    // Keep one independent feature per element so two MASS elements attached to
    // the same node remain two additive concentrated masses.
    auto region = std::make_shared<model::NodeRegion>(
        "__ABQ_MASS_" + std::to_string(element_id));
    region->add(point->nodes()[0]);

    auto point_mass = std::make_shared<feature::PointMass>();
    point_mass->region_ = std::move(region);
    point_mass->mass_   = mass;
    data.features.push_back(std::move(point_mass));
}

std::pair<Precision, std::string> ParserAbq::resolve_load_amplitude(const std::string& amplitude) {
    auto& data = *model()._data;
    auto* loadcase = active_loadcase();
    logging::error(loadcase != nullptr,
        "Cannot resolve a load amplitude without an active load case");
    const std::string procedure = loadcase->type_name();

    if (!amplitude.empty()) {
        logging::error(data.amplitudes.has(amplitude),
            "Unknown Abaqus amplitude '", amplitude, "'");

        if (procedure == "LINEARTRANSIENT" || procedure == "LINEARHARMONIC") {
            return {Precision(1), amplitude};
        }
        if (procedure == "LINEARSTATIC" || procedure == "LINEARBUCKLING") {
            return {data.amplitudes.get(amplitude)->evaluate(m_abq_state.step_period), std::string{}};
        }
        logging::error(procedure != "NONLINEARSTATIC",
            "Named load AMPLITUDE is not supported for nonlinear static/Riks proportional loading");
    }

    if (procedure == "LINEARTRANSIENT" && m_abq_state.step_amplitude == "RAMP") {
        const std::string name = "__ABQ_STEP_DEFAULT_AMPLITUDE";
        if (!data.amplitudes.has(name)) {
            auto generated = std::make_shared<bc::Amplitude>(name, bc::Interpolation::Linear);
            generated->add_sample(Precision(0), Precision(0));
            generated->add_sample(m_abq_state.step_period, Precision(1));
            model().add_amplitude(std::move(generated));
        }
        return {Precision(1), name};
    }

    logging::error(!(procedure == "NONLINEARSTATIC" && m_abq_state.step_amplitude == "STEP"),
        "STEP, AMPLITUDE=STEP cannot be represented by FEMaster nonlinear proportional load control");
    return {Precision(1), std::string{}};
}

void ParserAbq::register_common_commands(io::dsl::Registry& registry,
                                         std::shared_ptr<bool> assembly_scope) {
    // Register the complete Abaqus subset once per parser pass. ActiveMode below
    // controls execution; keeping one common grammar preserves scope handling
    // even when a command is only consumed in the current pass.
    commands::register_heading(registry);
    commands::register_part(registry, model());
    commands::register_end_part(registry, model());
    commands_abq::register_assembly(registry, assembly_scope);
    commands_abq::register_end_assembly(registry, assembly_scope);
    commands_abq::register_instance(registry, model());
    commands_abq::register_end_instance(registry);
    commands::register_nset(registry, model(), assembly_scope);
    commands::register_elset(registry, model(), assembly_scope);
    commands::register_surface(registry, model(), assembly_scope);
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
    commands_abq::register_mass(registry, *this, assembly_scope);
    commands_abq::register_equation(registry, *this);
    commands_abq::register_step(registry, *this);
    commands_abq::register_cload(registry, *this);
    commands_abq::register_boundary(registry, *this);
    commands_abq::register_dload(registry, *this);
    commands_abq::register_dsload(registry, *this);
}

void ParserAbq::configure_definition_pass(io::dsl::Registry& registry) {
    // Definitions may be referenced by part sections during the topology pass.
    // They therefore execute before any part or instance construction.
    auto assembly_scope = std::make_shared<bool>(false);
    commands::register_node(registry, model(), assembly_scope);
    commands_abq::register_element(registry, model(), assembly_scope);
    register_common_commands(registry, assembly_scope);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    for (const char* command : {
        "MATERIAL", "DENSITY", "ELASTIC", "HYPERELASTIC", "EXPANSION", "ORIENTATION"
    }) {
        registry.set_active_mode(command, io::dsl::ActiveMode::Active);
    }
}

void ParserAbq::configure_topology_pass(io::dsl::Registry& registry) {
    // The topology pass builds semantic parts, instances and part-local MASS
    // features. Assembly-level MASS records are consumed until compiled ELSETs
    // become available in the following pass.
    m_abq_state = ParserAbqState{};
    auto assembly_scope = std::make_shared<bool>(false);
    commands::register_node(registry, model(), assembly_scope);
    commands_abq::register_element(registry, model(), assembly_scope);
    register_common_commands(registry, assembly_scope);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    for (const char* command : {
        "PART", "ENDPART", "ASSEMBLY", "ENDASSEMBLY", "INSTANCE", "ENDINSTANCE",
        "NODE", "ELEMENT", "NSET", "ELSET", "SURFACE", "MASS", "AMPLITUDE", "SOLIDSECTION", "SHELLSECTION"
    }) {
        registry.set_active_mode(command, io::dsl::ActiveMode::Active);
    }
}

void ParserAbq::configure_assembly_pass(io::dsl::Registry& registry) {
    // Part-local point masses were already expanded by Model::compile(). This
    // pass only materializes assembly-local regions, MASS properties and nodal
    // transforms against dense identifiers.
    auto assembly_scope = std::make_shared<bool>(false);
    commands::register_node(registry, model(), assembly_scope);
    commands_abq::register_element(registry, model(), assembly_scope);
    register_common_commands(registry, assembly_scope);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    for (const char* command : {
        "ASSEMBLY", "ENDASSEMBLY", "NSET", "ELSET", "SURFACE", "MASS", "TRANSFORM"
    }) {
        registry.set_active_mode(command, io::dsl::ActiveMode::Active);
    }
}

void ParserAbq::configure_analysis_pass(io::dsl::Registry& registry) {
    // Analysis procedures and boundary/load data consume the final compiled
    // assembly, including all post-compile sets and nodal transforms.
    auto assembly_scope = std::make_shared<bool>(false);
    commands::register_node(registry, model(), assembly_scope);
    commands_abq::register_element(registry, model(), assembly_scope);
    register_common_commands(registry, assembly_scope);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    for (const char* command : {
        "EQUATION", "STEP", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS",
        "CLOAD", "BOUNDARY", "DLOAD", "DSLOAD", "ENDSTEP"
    }) {
        registry.set_active_mode(command, io::dsl::ActiveMode::Active);
    }
}

} // namespace fem::io::reader
