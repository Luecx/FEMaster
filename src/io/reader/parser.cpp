/**
 * @file parser.cpp
 * @brief Implements dependency-ordered FEMaster input-deck parsing.
 *
 * A complete command grammar is registered independently for every pass. The
 * registry then consumes inactive commands to preserve deck nesting while only
 * executing callbacks whose model dependencies have already been established.
 *
 * The implementation deliberately keeps the four passes explicit: definitions,
 * semantic topology, compiled assembly materialization and analysis. This makes
 * the one-way `Model::compile()` boundary and the state available to every
 * command visible in the top-level parser lifecycle.
 *
 * @see Parser
 * @see io::dsl::Registry
 * @see io::dsl::Engine
 * @see model::Model::compile
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "parser.h"

#include "../../core/logging.h"
#include "../../loadcase/loadcase.h"
#include "../../model/model.h"
#include "../dsl/engine.h"
#include "../dsl/file.h"
#include "../writer/writers.h"
#include "commands/register_amplitude.inl"
#include "commands/register_assembly.inl"
#include "commands/register_beam_section.inl"
#include "commands/register_cload.inl"
#include "commands/register_connector.inl"
#include "commands/register_contact.inl"
#include "commands/register_coupling.inl"
#include "commands/register_density.inl"
#include "commands/register_dload.inl"
#include "commands/register_elastic.inl"
#include "commands/register_element.inl"
#include "commands/register_elset.inl"
#include "commands/register_end_assembly.inl"
#include "commands/register_end_instance.inl"
#include "commands/register_end_part.inl"
#include "commands/register_equation.inl"
#include "commands/register_field.inl"
#include "commands/register_heading.inl"
#include "commands/register_hyperelastic.inl"
#include "commands/register_inertialload.inl"
#include "commands/register_instance.inl"
#include "commands/register_loadcase_begin.inl"
#include "commands/register_loadcase_constraintmethod.inl"
#include "commands/register_loadcase_constraintsummary.inl"
#include "commands/register_loadcase_damping.inl"
#include "commands/register_loadcase_frequency.inl"
#include "commands/register_loadcase_inertiarelief.inl"
#include "commands/register_loadcase_initialvelocity.inl"
#include "commands/register_loadcase_loads.inl"
#include "commands/register_loadcase_newmark.inl"
#include "commands/register_loadcase_nonlinear.inl"
#include "commands/register_loadcase_numeigenvalues.inl"
#include "commands/register_loadcase_rebalance.inl"
#include "commands/register_loadcase_request_stgeom.inl"
#include "commands/register_loadcase_request_stiffness.inl"
#include "commands/register_loadcase_sigma.inl"
#include "commands/register_loadcase_solver.inl"
#include "commands/register_loadcase_supports.inl"
#include "commands/register_loadcase_time.inl"
#include "commands/register_loadcase_topodensity.inl"
#include "commands/register_loadcase_topoexponent.inl"
#include "commands/register_loadcase_topoorient.inl"
#include "commands/register_loadcase_write_every.inl"
#include "commands/register_mass.inl"
#include "commands/register_material.inl"
#include "commands/register_node.inl"
#include "commands/register_normal.inl"
#include "commands/register_nset.inl"
#include "commands/register_orientation.inl"
#include "commands/register_overview.inl"
#include "commands/register_part.inl"
#include "commands/register_pload.inl"
#include "commands/register_point_mass.inl"
#include "commands/register_profile.inl"
#include "commands/register_rbm.inl"
#include "commands/register_rotary_inertia.inl"
#include "commands/register_sfset.inl"
#include "commands/register_shell_section.inl"
#include "commands/register_solid_section.inl"
#include "commands/register_spring.inl"
#include "commands/register_support.inl"
#include "commands/register_surface.inl"
#include "commands/register_thermalexpansion.inl"
#include "commands/register_tie.inl"
#include "commands/register_tload.inl"
#include "commands/register_truss_section.inl"
#include "commands/register_vload.inl"

#include <iostream>
#include <memory>
#include <utility>

namespace fem::io::reader {

/**
 * Constructs an idle parser and prepares the complete documentation grammar.
 *
 * Parsing itself starts with a fresh Model in `run()`. The initial model exists
 * so command callbacks required for registry documentation can be registered
 * immediately without executing a deck.
 */
Parser::Parser()
    : model_(std::make_shared<model::Model>()),
      writer_("") {
    configure_documentation_registry();
}

Parser::~Parser() = default;

/**
 * Evaluates a FEMaster deck through four dependency-ordered complete-file passes.
 *
 * Every invocation resets model and load-case construction state. The definition
 * pass first creates resources referenced by sections. The topology pass then
 * builds sparse semantic Parts and Instances before the explicit one-way compile
 * operation creates dense assembly storage. The assembly pass materializes
 * compiled regions, fields and shell normals, and the analysis pass finally
 * executes loads, constraints and load cases while writing results.
 *
 * @param input_path Input deck evaluated by every parser pass.
 * @param output_path Optional base path for result files. The input path is used
 *                    when this value is empty.
 * @param writer_formats Result formats enabled for analysis output.
 */
void Parser::run(const std::string& input_path,
                 const std::string& output_path,
                 const io::writer::WriterFileFormats& writer_formats) {
    // Reset all mutable state so each run represents an independent deck
    model_ = std::make_shared<model::Model>();
    active_loadcase_.reset();
    next_loadcase_id_ = 1;

    // Collect shared definitions before sections resolve their dependencies
    run_definition_pass(input_path);

    // Construct semantic Parts, Instances and part-local topology
    run_topology_pass(input_path);

    // Cross the one-way boundary into deterministic dense assembly storage
    model_->compile();

    // Materialize post-compile assembly regions, fields and shell normals
    run_assembly_pass(input_path);

    // Execute loads, constraints and load cases against the completed assembly
    run_analysis_pass(input_path, output_path, writer_formats);
}

/**
 * Prints one requested view of the registered FEMaster command grammar.
 *
 * Documentation operates on the persistent registry prepared independently of
 * the temporary registries used for actual deck passes. Unsupported output
 * formats currently fall back to the text renderer without changing the
 * requested action.
 *
 * @param opts Documentation action, query and formatting selection.
 */
void Parser::document(const DocOptions& opts) const {
    using A = DocOptions::Action;
    using F = DocOptions::Format;
    using V = DocOptions::Verbosity;

    if (opts.format != F::Text) {
        std::cout << "(Note) Only text output is implemented currently. Falling back to text.\n\n";
    }

    // Dispatch the selected documentation view to the complete grammar registry
    switch (opts.action) {
        case A::List:       documentation_registry_.print_index(); break;
        case A::Show:       documentation_registry_.print_help(opts.cmd, opts.verbosity == V::Compact); break;
        case A::Tokens:     documentation_registry_.print_tokens(opts.cmd); break;
        case A::Variants:   documentation_registry_.print_variants(opts.cmd); break;
        case A::Search:     documentation_registry_.print_search(opts.query, opts.regex); break;
        case A::WhereToken: documentation_registry_.print_where_token(opts.query); break;
        case A::All:        documentation_registry_.print_help({}, false); break;
    }
}

/**
 * Returns the mutable model currently constructed or analyzed by the parser.
 *
 * The model is replaced at the beginning of every `run()` invocation. Access
 * is rejected when no current model exists so command callbacks never operate
 * on an invalid parser context.
 *
 * @return Mutable current model.
 */
model::Model& Parser::model() {
    logging::error(model_ != nullptr,
        "Parser: model is not initialized");
    return *model_;
}

/**
 * Returns the model currently constructed or analyzed by the parser.
 *
 * The same initialization invariant as the mutable overload is enforced for
 * read-only parser and documentation operations.
 *
 * @return Read-only current model.
 */
const model::Model& Parser::model() const {
    logging::error(model_ != nullptr,
        "Parser: model is not initialized");
    return *model_;
}

const io::dsl::Registry& Parser::registry() const { return documentation_registry_; }

/**
 * Activates a load case and supplies its parser-owned analysis context.
 *
 * Activation assigns the next sequential load-case identifier together with
 * the current result writer and compiled model. The parser then retains sole
 * ownership while consecutive load-case commands configure the analysis.
 *
 * @param loadcase Newly constructed concrete load case.
 */
void Parser::begin_loadcase(loadcase::LoadCase::Ptr loadcase) {
    logging::error(active_loadcase_ == nullptr,
        "Parser: nested load cases are not supported");
    logging::error(loadcase != nullptr,
        "Parser: cannot activate a null load case");
    logging::error(model_ != nullptr,
        "Parser: cannot activate a load case without a model");

    loadcase->id     = next_loadcase_id_++;
    loadcase->writer = &writer_;
    loadcase->model  = model_.get();
    active_loadcase_ = std::move(loadcase);
}

/**
 * Completes and executes the active load case.
 *
 * Ownership is removed from the parser before execution so the definition
 * scope is already closed while the solver runs. The load case is destroyed
 * automatically after the analysis finishes or propagates an exception.
 */
void Parser::end_loadcase() {
    logging::error(active_loadcase_ != nullptr,
        "Parser: cannot end a load case when none is active");

    auto loadcase = std::move(active_loadcase_);
    loadcase->run();
}

/**
 * Returns the load case currently configured by consecutive parser commands.
 *
 * The returned pointer is non-owning and remains valid only until
 * `end_loadcase()` executes and releases the active analysis.
 *
 * @return Active load case, or `nullptr` outside a load-case scope.
 */
loadcase::LoadCase* Parser::active_loadcase() {
    return active_loadcase_.get();
}

/**
 * Executes the definition pass over the complete input deck.
 *
 * A fresh registry recognizes the entire FEMaster grammar but activates only
 * shared definitions required by later topology and section commands.
 *
 * @param input_path Input deck replayed for definition callbacks.
 */
void Parser::run_definition_pass(const std::string& input_path) {
    // Configure an independent grammar and replay the complete deck once
    io::dsl::Registry registry;
    configure_definition_pass(registry);
    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);
}

/**
 * Executes the semantic-topology pass over the complete input deck.
 *
 * Parts, Instances, local topology and section assignments are constructed while
 * dense assembly-dependent callbacks remain consume-only.
 *
 * @param input_path Input deck replayed for topology callbacks.
 */
void Parser::run_topology_pass(const std::string& input_path) {
    // Configure an independent grammar and replay the complete deck once
    io::dsl::Registry registry;
    configure_topology_pass(registry);
    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);
}

/**
 * Executes post-compile assembly materialization over the complete input deck.
 *
 * Assembly regions, fields and prescribed normals can now resolve dense global
 * identifiers and compiled element-domain dimensions. After parsing, missing
 * shell element-nodal normals are reconstructed and equalized so all subsequent
 * analyses see a complete reference-normal field.
 *
 * @param input_path Input deck replayed for assembly callbacks.
 */
void Parser::run_assembly_pass(const std::string& input_path) {
    // Materialize commands whose storage or identifiers require compilation
    io::dsl::Registry registry;
    configure_assembly_pass(registry);
    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);

    // Complete solver-facing shell reference data before any load case executes
    model_->build_shell_element_normals();
}

/**
 * Executes loads, constraints and load cases over the completed assembly.
 *
 * Result writers are initialized from the requested output base and receive the
 * compiled model before analysis callbacks begin. A final complete-deck pass then
 * executes only commands not consumed by the preceding construction passes. The
 * writers are closed after the last load case finishes.
 *
 * @param input_path Input deck replayed for analysis callbacks and used as the
 *                   default result base.
 * @param output_path Optional explicit result base or result filename.
 * @param writer_formats Result formats enabled for analysis output.
 */
void Parser::run_analysis_pass(const std::string& input_path,
                               const std::string& output_path,
                               const io::writer::WriterFileFormats& writer_formats) {
    // Derive a format-independent writer base by removing a known result suffix
    std::string writer_base = output_path.empty() ? input_path : output_path;
    for (const std::string& ext : {std::string(".res"), std::string(".frd"), std::string(".inp")}) {
        if (writer_base.size() >= ext.size()
         && writer_base.compare(writer_base.size() - ext.size(), ext.size(), ext) == 0) {
            writer_base.resize(writer_base.size() - ext.size());
            break;
        }
    }

    // Publish the completed model before load-case commands start producing frames
    writer_ = io::writer::ResultWriters(writer_base, writer_formats);
    writer_.write_model_data(*model_->_data);

    // Execute the final pass against the fully materialized assembly
    io::dsl::Registry registry;
    configure_analysis_pass(registry);
    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);

    // Flush and close every enabled result format after the final command
    writer_.close();
}

/**
 * Configures the pass that collects shared model definitions.
 *
 * The complete grammar remains available for nesting and validation, but only
 * material laws, profiles and coordinate orientations execute callbacks. These
 * definitions may therefore be referenced later regardless of deck order.
 *
 * @param registry Empty pass-local registry to configure.
 */
void Parser::configure_definition_pass(io::dsl::Registry& registry) {
    // Register the complete grammar while suppressing all callbacks by default
    register_commands(registry);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);

    // Activate global definitions required by later section construction
    registry.set_active_mode("MATERIAL", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELASTIC", io::dsl::ActiveMode::Active);
    registry.set_active_mode("HYPERELASTIC", io::dsl::ActiveMode::Active);
    registry.set_active_mode("DENSITY", io::dsl::ActiveMode::Active);
    registry.set_active_mode("THERMALEXPANSION", io::dsl::ActiveMode::Active);
    registry.set_active_mode("PROFILE", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ORIENTATION", io::dsl::ActiveMode::Active);
}

/**
 * Configures construction of sparse semantic topology before compilation.
 *
 * Part and Instance scopes, local entities, regions and part-level section
 * assignments execute while global assembly fields and analyses are consumed
 * without mutating state.
 *
 * @param registry Empty pass-local registry to configure.
 */
void Parser::configure_topology_pass(io::dsl::Registry& registry) {
    // Register the complete grammar while suppressing all callbacks by default
    register_commands(registry);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);

    // Activate commands that construct semantic Parts and Instances
    for (const char* command : {
        "PART", "ENDPART", "ASSEMBLY", "ENDASSEMBLY", "INSTANCE", "ENDINSTANCE",
        "NODE", "ELEMENT", "NSET", "ELSET", "SURFACE", "SFSET",
        "SOLIDSECTION", "BEAMSECTION", "TRUSSSECTION", "SHELLSECTION",
        "MASS", "ROTARYINERTIA", "SPRING",
    }) {
        registry.set_active_mode(command, io::dsl::ActiveMode::Active);
    }
}

/**
 * Configures materialization of assembly entities after model compilation.
 *
 * Assembly scope is replayed so global regions and surfaces can resolve Instance
 * maps. Fields and prescribed normals may now use compiled domain dimensions and
 * dense global identifiers. Point-element properties are also materialized here
 * so assembly-level ELSET assignments exist before result writers are initialized.
 *
 * @param registry Empty pass-local registry to configure.
 */
void Parser::configure_assembly_pass(io::dsl::Registry& registry) {
    // Register the complete grammar while suppressing all callbacks by default
    register_commands(registry);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);

    // Activate post-compile assembly, point-property and field materialization commands
    for (const char* command : {
        "ASSEMBLY", "ENDASSEMBLY", "NSET", "ELSET", "SURFACE", "SFSET", "FIELD", "NORMAL",
        "MASS", "ROTARYINERTIA", "SPRING"
    }) {
        registry.set_active_mode(command, io::dsl::ActiveMode::Active);
    }
}

/**
 * Configures the final load, constraint and load-case execution pass.
 *
 * Commands are active by default because this pass owns all remaining model and
 * analysis operations. Definitions, topology and assembly materialization are
 * switched to consume-only so their callbacks cannot duplicate established
 * state while their scopes remain syntactically available.
 *
 * @param registry Empty pass-local registry to configure.
 */
void Parser::configure_analysis_pass(io::dsl::Registry& registry) {
    // Register the complete grammar and activate remaining callbacks by default
    register_commands(registry);
    registry.set_active_mode(io::dsl::ActiveMode::Active);

    // Consume commands whose state was finalized by an earlier parser pass
    for (const char* command : {
        "PART", "ENDPART", "ASSEMBLY", "ENDASSEMBLY", "INSTANCE", "ENDINSTANCE",
        "NODE", "ELEMENT", "NSET", "ELSET", "SURFACE", "SFSET",
        "MATERIAL", "ELASTIC", "HYPERELASTIC", "DENSITY", "THERMALEXPANSION",
        "PROFILE", "ORIENTATION", "SOLIDSECTION", "BEAMSECTION", "TRUSSSECTION",
        "SHELLSECTION", "FIELD", "NORMAL", "MASS", "ROTARYINERTIA", "SPRING"
    }) {
        registry.set_active_mode(command, io::dsl::ActiveMode::ConsumeOnly);
    }
}

/**
 * Rebuilds the persistent registry used for command-language documentation.
 *
 * All registered commands are active so index, token, variant and search output
 * can inspect the complete FEMaster grammar independently of parser-pass masks.
 */
void Parser::configure_documentation_registry() {
    // Replace any previous grammar and expose every command to documentation
    documentation_registry_ = io::dsl::Registry{};
    register_commands(documentation_registry_);
    documentation_registry_.set_active_mode(io::dsl::ActiveMode::Active);
}

/**
 * Registers the complete FEMaster command grammar for one independent registry.
 *
 * Every parsing pass uses the same grammar. Active modes control whether a
 * command executes its callback or is only consumed to preserve the original
 * scope hierarchy. Parent-dependent behavior is represented directly by DSL
 * conditions and variants rather than duplicated parser-local scope flags.
 *
 * Model-oriented commands receive the current Model directly. Load-case commands
 * receive the Parser because they coordinate consecutive commands through the
 * active load-case state and result writers.
 *
 * @param registry Registry that receives the complete command definitions.
 */
void Parser::register_commands(io::dsl::Registry& registry) {
    // Validate the callback target before registering model-oriented commands
    logging::error(model_ != nullptr,
        "Parser: model must exist before registering commands");

    auto& mdl = *model_;

    // Register semantic Part/Instance topology and scope-aware assembly commands
    commands::register_part        (registry, mdl);
    commands::register_end_part    (registry, mdl);
    commands::register_assembly    (registry);
    commands::register_end_assembly(registry);
    commands::register_instance    (registry, mdl);
    commands::register_end_instance(registry);
    commands::register_node        (registry, mdl);
    commands::register_element     (registry, mdl);
    commands::register_nset        (registry, mdl);
    commands::register_elset       (registry, mdl);
    commands::register_surface     (registry, mdl);
    commands::register_sfset       (registry, mdl);

    // Register the root model and load-case terminator scopes
    registry.command("MODEL", [](io::dsl::Command& command) {
        command.allow_if(io::dsl::Condition::parent_is("ROOT"));
        command.keyword(io::dsl::KeywordSpec::make().key("NAME").optional());
        command.variant(io::dsl::Variant::make());
    });
    registry.command("END", [](io::dsl::Command& command) {
        command.allow_if(io::dsl::Condition::parent_is("LOADCASE"));
        command.closes_parent();
        command.variant(io::dsl::Variant::make());
    });

    // Register field, material, profile and section definitions
    commands::register_heading(registry);
    commands::register_field(registry, mdl);
    commands::register_normal(registry, mdl);
    commands::register_material(registry, mdl);
    commands::register_elastic(registry, mdl);
    commands::register_hyperelastic(registry, mdl);
    commands::register_density(registry, mdl);
    commands::register_thermal_expansion(registry, mdl);
    commands::register_orientation(registry, mdl);
    commands::register_profile(registry, mdl);
    commands::register_solid_section(registry, mdl);
    commands::register_beam_section(registry, mdl);
    commands::register_truss_section(registry, mdl);
    commands::register_shell_section(registry, mdl);
    commands::register_mass(registry, mdl);
    commands::register_rotary_inertia(registry, mdl);
    commands::register_spring(registry, mdl);

    // Register loads, constraints, features and model diagnostics
    commands::register_cload(registry, mdl);
    commands::register_dload(registry, mdl);
    commands::register_pload(registry, mdl);
    commands::register_tload(registry, mdl);
    commands::register_vload(registry, mdl);
    commands::register_inertialload(registry, mdl);
    commands::register_rbm(registry, mdl);
    commands::register_support(registry, mdl);
    commands::register_amplitude(registry, mdl);
    commands::register_connector(registry, mdl);
    commands::register_coupling(registry, mdl);
    commands::register_tie(registry, mdl);
    commands::register_contact(registry, mdl);
    commands::register_point_mass(registry, mdl);
    commands::register_overview(registry, mdl);
    commands::register_equation(registry, mdl);

    // Register load-case creation, solver settings and result requests
    commands::register_loadcase_begin(registry, *this);
    commands::register_loadcase_supports(registry, *this);
    commands::register_loadcase_loads(registry, *this);
    commands::register_loadcase_solver(registry, *this);
    commands::register_loadcase_constraintmethod(registry, *this);
    commands::register_loadcase_frequency(registry, *this);
    commands::register_loadcase_request_stiffness(registry, *this);
    commands::register_loadcase_request_stgeom(registry, *this);
    commands::register_loadcase_numeigenvalues(registry, *this);
    commands::register_loadcase_sigma(registry, *this);
    commands::register_loadcase_topodensity(registry, *this);
    commands::register_loadcase_topoorient(registry, *this);
    commands::register_loadcase_topoexponent(registry, *this);
    commands::register_loadcase_constraintsummary(registry, *this);
    commands::register_loadcase_nonlinear(registry, *this);
    commands::register_loadcase_time(registry, *this);
    commands::register_loadcase_write_every(registry, *this);
    commands::register_loadcase_damping(registry, *this);
    commands::register_loadcase_newmark(registry, *this);
    commands::register_loadcase_initialvelocity(registry, *this);
    commands::register_loadcase_inertiarelief(registry, *this);
    commands::register_loadcase_rebalance(registry, *this);
}

} // namespace fem::io::reader
