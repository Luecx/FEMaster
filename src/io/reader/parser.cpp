/**
 * @file parser.cpp
 * @brief Implements one-shot parsing and explicit FEMaster deck processing.
 *
 * The input file is parsed exactly once into `dsl::Deck`. Semantic model
 * construction then follows a visible dependency order: global definitions,
 * Part-local topology, Instances, `Model::compile()`, assembly materialization,
 * compiled model features and finally load cases. Scope-sensitive commands are
 * selected from their stored parent nodes instead of being replayed in parser
 * stages.
 *
 * The registered `.inl` command callbacks remain the semantic implementation.
 * `ParsedCommand::enter()`, `execute()` and `leave()` control when those existing
 * callbacks run; syntax validation, variant selection and data aggregation have
 * already completed in `dsl::DeckParser`.
 *
 * @see Parser
 * @see io::dsl::Deck
 * @see io::dsl::DeckParser
 * @see model::Model::compile
 *
 * @author Finn Eggers
 * @date 26.08.2026
 */

#include "parser.h"

#include "../../core/logging.h"
#include "../../loadcase/loadcase.h"
#include "../../model/model.h"
#include "../dsl/deck_parser.h"
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
 * Constructs an idle parser and prepares the native FEMaster documentation grammar.
 */
Parser::Parser()
    : model_(std::make_shared<model::Model>()),
      writer_("") {
    configure_documentation_registry();
}

Parser::~Parser() = default;

/**
 * Parses one complete input deck and executes its semantic commands afterwards.
 *
 * Every run owns one registry for the entire parse/process lifetime. This is
 * required because parsed command and segment pointers refer directly to that
 * registry. The file itself is consumed only once; all later dependency ordering
 * operates on the in-memory deck representation.
 *
 * @param input_path Input deck parsed exactly once.
 * @param output_path Optional base path for result files.
 * @param writer_formats Result formats enabled for analysis output.
 */
void Parser::run(const std::string& input_path,
                 const std::string& output_path,
                 const io::writer::WriterFileFormats& writer_formats) {
    // Reset all mutable state so each run represents an independent deck.
    model_ = std::make_shared<model::Model>();
    active_loadcase_.reset();
    next_loadcase_id_ = 1;

    // Register the complete grammar once and parse the complete source once.
    io::dsl::Registry registry;
    register_commands(registry);

    io::dsl::File       file(input_path);
    io::dsl::DeckParser deck_parser(registry);
    const io::dsl::Deck deck = deck_parser.parse(file);

    // Semantic dependencies are now explicit and independent of source order.
    process_deck(deck, input_path, output_path, writer_formats);
}

/**
 * Executes the parsed native FEMaster deck in explicit model-dependency order.
 *
 * The function deliberately spells out the semantic order from top to bottom.
 * Loops only iterate repeated input scopes such as MATERIAL, PART, ASSEMBLY,
 * COUPLING and LOADCASE; command ordering inside each scope remains visible at
 * the call site. `Model::compile()` forms the one-way boundary between sparse
 * Part/Instance topology and dense assembly materialization.
 */
void Parser::process_deck(const io::dsl::Deck&                  deck,
                          const std::string&                    input_path,
                          const std::string&                    output_path,
                          const io::writer::WriterFileFormats& writer_formats) {
    const auto& root = deck.root();

    // ---------------------------------------------------------------------
    // Global definitions required by sections and later load definitions
    // ---------------------------------------------------------------------
    root.execute_children("HEADING");

    for (const auto* material : root.children("MATERIAL")) {
        material->enter();

        material->execute_children("ELASTIC");
        material->execute_children("HYPERELASTIC");
        material->execute_children("DENSITY");
        material->execute_children("THERMALEXPANSION");

        material->leave();
    }

    root.execute_children("PROFILE");
    root.execute_children("ORIENTATION");
    root.execute_children("AMPLITUDE");

    // ---------------------------------------------------------------------
    // Default-Part topology before Model::compile()
    // ---------------------------------------------------------------------
    root.execute_children("NODE");
    root.execute_children("ELEMENT");

    root.execute_children("NSET");
    root.execute_children("ELSET");
    root.execute_children("SURFACE");
    root.execute_children("SFSET");

    root.execute_children("SOLIDSECTION");
    root.execute_children("BEAMSECTION");
    root.execute_children("TRUSSSECTION");
    root.execute_children("SHELLSECTION");

    root.execute_children("MASS");
    root.execute_children("ROTARYINERTIA");
    root.execute_children("SPRING");

    // ---------------------------------------------------------------------
    // Explicit Part topology before Model::compile()
    // ---------------------------------------------------------------------
    for (const auto* part : root.children("PART")) {
        part->enter();

        part->execute_children("NODE");
        part->execute_children("ELEMENT");

        part->execute_children("NSET");
        part->execute_children("ELSET");
        part->execute_children("SURFACE");
        part->execute_children("SFSET");

        part->execute_children("SOLIDSECTION");
        part->execute_children("BEAMSECTION");
        part->execute_children("TRUSSSECTION");
        part->execute_children("SHELLSECTION");

        part->execute_children("MASS");
        part->execute_children("ROTARYINERTIA");
        part->execute_children("SPRING");

        part->leave();
    }

    // ---------------------------------------------------------------------
    // Instances depend on completed Parts
    // ---------------------------------------------------------------------
    root.execute_children("INSTANCE");

    for (const auto* assembly : root.children("ASSEMBLY")) {
        assembly->execute_children("INSTANCE");
    }

    // ---------------------------------------------------------------------
    // Transition from sparse topology to dense assembly data
    // ---------------------------------------------------------------------
    model_->compile();

    // ---------------------------------------------------------------------
    // Assembly regions, properties and dense fields
    // ---------------------------------------------------------------------
    root.execute_children("FIELD");
    root.execute_children("NORMAL");

    for (const auto* assembly : root.children("ASSEMBLY")) {
        assembly->execute_children("NSET");
        assembly->execute_children("ELSET");
        assembly->execute_children("SURFACE");
        assembly->execute_children("SFSET");

        assembly->execute_children("MASS");
        assembly->execute_children("ROTARYINERTIA");
        assembly->execute_children("SPRING");

        assembly->execute_children("FIELD");
        assembly->execute_children("NORMAL");
    }

    // Complete reference normals before constraints or analyses consume shells.
    model_->build_shell_element_normals();

    // ---------------------------------------------------------------------
    // Compiled model features
    // ---------------------------------------------------------------------
    root.execute_children("POINTMASS");

    // ---------------------------------------------------------------------
    // Root-level collectors and constraints
    // ---------------------------------------------------------------------
    root.execute_children("SUPPORT");

    root.execute_children("CLOAD");
    root.execute_children("DLOAD");
    root.execute_children("PLOAD");
    root.execute_children("TLOAD");
    root.execute_children("VLOAD");
    root.execute_children("INERTIALLOAD");

    root.execute_children("RBM");
    root.execute_children("CONNECTOR");
    root.execute_children("TIE");
    root.execute_children("CONTACT");
    root.execute_children("EQUATION");

    // ---------------------------------------------------------------------
    // Assembly-level collectors and constraints
    // ---------------------------------------------------------------------
    for (const auto* assembly : root.children("ASSEMBLY")) {
        assembly->execute_children("SUPPORT");

        assembly->execute_children("CLOAD");
        assembly->execute_children("DLOAD");
        assembly->execute_children("PLOAD");
        assembly->execute_children("TLOAD");
        assembly->execute_children("VLOAD");
        assembly->execute_children("INERTIALLOAD");

        assembly->execute_children("RBM");
        assembly->execute_children("CONNECTOR");
        assembly->execute_children("TIE");
        assembly->execute_children("CONTACT");
        assembly->execute_children("EQUATION");
    }

    // ---------------------------------------------------------------------
    // Couplings
    // ---------------------------------------------------------------------
    for (const auto* coupling : root.children("COUPLING")) {
        coupling->enter();

        coupling->execute_children("KINEMATIC");
        coupling->execute_children("DISTRIBUTING");

        coupling->leave();
    }

    for (const auto* assembly : root.children("ASSEMBLY")) {
        for (const auto* coupling : assembly->children("COUPLING")) {
            coupling->enter();

            coupling->execute_children("KINEMATIC");
            coupling->execute_children("DISTRIBUTING");

            coupling->leave();
        }
    }

    root.execute_children("OVERVIEW");

    // ---------------------------------------------------------------------
    // Result writers and load-case execution
    // ---------------------------------------------------------------------
    initialize_writers(input_path, output_path, writer_formats);

    for (const auto* loadcase : root.children("LOADCASE")) {
        loadcase->enter();
        loadcase->execute_children();
        loadcase->leave();
    }

    close_writers();
}

/**
 * Initializes enabled result writers from the requested output path and publishes
 * the completely materialized model before any load case starts writing frames.
 */
void Parser::initialize_writers(const std::string&                    input_path,
                                const std::string&                    output_path,
                                const io::writer::WriterFileFormats& writer_formats) {
    std::string writer_base = output_path.empty() ? input_path : output_path;
    for (const std::string& ext : {std::string(".res"), std::string(".frd"), std::string(".inp")}) {
        if (writer_base.size() >= ext.size()
         && writer_base.compare(writer_base.size() - ext.size(), ext.size(), ext) == 0) {
            writer_base.resize(writer_base.size() - ext.size());
            break;
        }
    }

    writer_ = io::writer::ResultWriters(writer_base, writer_formats);
    writer_.write_model_data(*model_->_data);
}

/**
 * Flushes and closes every enabled result writer after semantic analysis processing.
 */
void Parser::close_writers() {
    writer_.close();
}

/**
 * Prints one requested view of the registered FEMaster command grammar.
 */
void Parser::document(const DocOptions& opts) const {
    using A = DocOptions::Action;
    using F = DocOptions::Format;
    using V = DocOptions::Verbosity;

    if (opts.format != F::Text) {
        std::cout << "(Note) Only text output is implemented currently. Falling back to text.\n\n";
    }

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
 */
model::Model& Parser::model() {
    logging::error(model_ != nullptr,
        "Parser: model is not initialized");
    return *model_;
}

/**
 * Returns the model currently constructed or analyzed by the parser.
 */
const model::Model& Parser::model() const {
    logging::error(model_ != nullptr,
        "Parser: model is not initialized");
    return *model_;
}

const io::dsl::Registry& Parser::registry() const { return documentation_registry_; }

/**
 * Activates a load case and supplies its parser-owned analysis context.
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
 */
void Parser::end_loadcase() {
    logging::error(active_loadcase_ != nullptr,
        "Parser: cannot end a load case when none is active");

    auto loadcase = std::move(active_loadcase_);
    loadcase->run();
}

/**
 * Returns the load case currently configured by consecutive parser commands.
 */
loadcase::LoadCase* Parser::active_loadcase() {
    return active_loadcase_.get();
}

/**
 * Rebuilds the persistent registry used for command-language documentation.
 */
void Parser::configure_documentation_registry() {
    documentation_registry_ = io::dsl::Registry{};
    register_commands(documentation_registry_);
}

/**
 * Registers the complete native FEMaster command grammar exactly once per run.
 *
 * Registration describes syntax and semantic callbacks only. No parser stage or
 * command activation mode is attached to the grammar; execution timing belongs
 * exclusively to `process_deck()`.
 */
void Parser::register_commands(io::dsl::Registry& registry) {
    logging::error(model_ != nullptr,
        "Parser: model must exist before registering commands");

    auto& mdl = *model_;

    // Semantic Part/Instance topology and scope-aware assembly commands
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

    // Root model and native load-case terminator scopes
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

    // Field, material, profile and section definitions
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

    // Loads, constraints, features and model diagnostics
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

    // Load-case creation, solver settings and result requests
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
