/**
 * @file parser_abq.cpp
 * @brief Implements one-shot Abaqus parsing and explicit semantic processing.
 *
 * The Abaqus grammar is registered once and the source file is parsed once by
 * the common `Parser::run()` lifecycle. This file then exposes the semantic
 * dependency order directly: definitions, Part topology, Instances,
 * `Model::compile()`, assembly regions/properties/transforms, constraints and
 * finally the Abaqus analysis step.
 *
 * Scope-sensitive records are selected from the stored parsed tree. No parser
 * stages or command activation modes are required, and syntax/variant matching
 * is never repeated after the initial parse.
 *
 * @author Finn Eggers
 * @date 26.08.2026
 */

#include "parser_abq.h"

#include "../../bc/amplitude.h"
#include "../../loadcase/loadcase.h"
#include "../../model/model.h"
#include "../dsl/registry.h"
#include "commands/register_assembly.inl"
#include "commands/register_density.inl"
#include "commands/register_elastic.inl"
#include "commands/register_elset.inl"
#include "commands/register_end_assembly.inl"
#include "commands/register_end_instance.inl"
#include "commands/register_end_part.inl"
#include "commands/register_equation.inl"
#include "commands/register_coupling.inl"
#include "commands/register_heading.inl"
#include "commands/register_hyperelastic.inl"
#include "commands/register_instance.inl"
#include "commands/register_mass.inl"
#include "commands/register_material.inl"
#include "commands/register_node.inl"
#include "commands/register_nset.inl"
#include "commands/register_part.inl"
#include "commands/register_rotary_inertia.inl"
#include "commands/register_spring.inl"
#include "commands/register_surface.inl"
#include "commands_abq/register_amplitude.inl"
#include "commands_abq/register_boundary.inl"
#include "commands_abq/register_cload.inl"
#include "commands_abq/register_dload.inl"
#include "commands_abq/register_dsload.inl"
#include "commands_abq/register_element.inl"
#include "commands_abq/register_expansion.inl"
#include "commands_abq/register_orientation.inl"
#include "commands_abq/register_shell_section.inl"
#include "commands_abq/register_solid_section.inl"
#include "commands_abq/register_step.inl"
#include "commands_abq/register_transform.inl"

#include <memory>
#include <string>
#include <utility>

namespace fem::io::reader {

ParserAbq::ParserAbq() {
    // The base constructor can only dispatch the native virtual implementation.
    // Rebuild documentation here after the complete Abaqus dynamic type exists.
    configure_documentation_registry();
}

ParserAbqState& ParserAbq::abaqus_state() { return m_abq_state; }
const ParserAbqState& ParserAbq::abaqus_state() const { return m_abq_state; }

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

/**
 * Registers command definitions shared by the supported Abaqus subset.
 */
void ParserAbq::register_common_commands(io::dsl::Registry& registry) {
    commands::register_heading(registry);
    commands::register_part(registry, model());
    commands::register_end_part(registry, model());
    commands::register_assembly(registry);
    commands::register_end_assembly(registry);
    commands::register_instance(registry, model());
    commands::register_end_instance(registry);
    commands::register_nset(registry, model());
    commands::register_elset(registry, model());
    commands::register_surface(registry, model());
    commands::register_material(registry, model());
    commands::register_density(registry, model());
    commands::register_elastic(registry, model());
    commands::register_hyperelastic(registry, model());
    commands::register_mass(registry, model());
    commands::register_rotary_inertia(registry, model());
    commands::register_spring(registry, model());
    commands::register_equation(registry, model());
    commands::register_coupling(registry, model());

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
}

/**
 * Registers the complete supported Abaqus grammar once for parsing and processing.
 */
void ParserAbq::register_commands(io::dsl::Registry& registry) {
    commands::register_node(registry, model());
    commands_abq::register_element(registry, model());
    register_common_commands(registry);
}

/**
 * Executes the parsed Abaqus deck in explicit model-dependency order.
 *
 * Part-local sections and point properties are assigned before compilation so
 * ordinary Instance expansion copies them into assembly space. Assembly NSET,
 * ELSET, SURFACE, point properties and TRANSFORM records execute only after the
 * dense mappings exist. The STEP subtree remains source-ordered because its
 * procedure card establishes one active load case consumed by following load and
 * boundary cards until ENDSTEP executes it.
 */
void ParserAbq::process_deck(const io::dsl::Deck&                  deck,
                             const std::string&                    input_path,
                             const std::string&                    output_path,
                             const io::writer::WriterFileFormats& writer_formats) {
    using NodeId = io::dsl::Deck::NodeId;

    m_abq_state = ParserAbqState{};

    // ---------------------------------------------------------------------
    // Global material resources, orientations and amplitudes
    // ---------------------------------------------------------------------
    deck.execute_roots("HEADING");

    for (const NodeId material : deck.roots("MATERIAL")) {
        deck.enter(material);
        deck.execute_children(material, {"ELASTIC", "HYPERELASTIC", "DENSITY", "EXPANSION"});
        deck.leave(material);
    }

    deck.execute_roots({"ORIENTATION", "AMPLITUDE"});

    // ---------------------------------------------------------------------
    // Default-Part and explicit Part topology before compilation
    // ---------------------------------------------------------------------
    const auto process_part_topology = [&deck](NodeId part) {
        deck.execute_children(part, {
            "NODE", "ELEMENT",
            "NSET", "ELSET", "SURFACE",
            "SOLIDSECTION", "SHELLSECTION",
            "MASS", "ROTARYINERTIA", "SPRING"
        });
    };

    deck.execute_roots({
        "NODE", "ELEMENT",
        "NSET", "ELSET", "SURFACE",
        "SOLIDSECTION", "SHELLSECTION",
        "MASS", "ROTARYINERTIA", "SPRING"
    });

    for (const NodeId part : deck.roots("PART")) {
        deck.enter(part);
        process_part_topology(part);
        deck.leave(part);
    }

    // Instances require complete Parts but are still semantic pre-compile objects.
    deck.execute_roots("INSTANCE");
    for (const NodeId assembly : deck.roots("ASSEMBLY")) {
        deck.execute_children(assembly, "INSTANCE");
    }

    // ---------------------------------------------------------------------
    // Flatten Parts and Instances into the dense assembly exactly once
    // ---------------------------------------------------------------------
    model().compile();

    // ---------------------------------------------------------------------
    // Assembly regions, point properties and nodal transforms
    // ---------------------------------------------------------------------
    deck.execute_roots("TRANSFORM");

    for (const NodeId assembly : deck.roots("ASSEMBLY")) {
        deck.execute_children(assembly, {
            "NSET", "ELSET", "SURFACE",
            "MASS", "ROTARYINERTIA", "SPRING", "TRANSFORM"
        });
    }

    model().build_shell_element_normals();

    // ---------------------------------------------------------------------
    // Compiled initial boundaries and global constraints
    // ---------------------------------------------------------------------
    deck.execute_roots("BOUNDARY");
    deck.execute_roots("EQUATION");
    for (const NodeId assembly : deck.roots("ASSEMBLY")) {
        deck.execute_children(assembly, {"BOUNDARY", "EQUATION"});
    }

    const auto process_coupling = [&deck](NodeId coupling) {
        deck.enter(coupling);
        deck.execute_children(coupling, {"KINEMATIC", "DISTRIBUTING"});
        deck.leave(coupling);
    };

    for (const NodeId coupling : deck.roots("COUPLING")) process_coupling(coupling);
    for (const NodeId assembly : deck.roots("ASSEMBLY")) {
        for (const NodeId coupling : deck.children(assembly, "COUPLING")) process_coupling(coupling);
    }

    // ---------------------------------------------------------------------
    // Result writers and final Abaqus STEP execution
    // ---------------------------------------------------------------------
    initialize_writers(input_path, output_path, writer_formats);

    for (const NodeId step : deck.roots("STEP")) {
        deck.enter(step);

        // Procedure, loads, boundaries, output requests and ENDSTEP retain their
        // original order because they operate on one active Abaqus step state.
        deck.execute_children(step);

        deck.leave(step);
    }

    close_writers();
}

} // namespace fem::io::reader
