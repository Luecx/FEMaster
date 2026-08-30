/**
 * @file parser.cpp
 * @brief Implements one-shot parsing and explicit FEMaster deck processing.
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "parser.h"

#include "../../core/logging.h"
#include "../../loadcase/loadcase.h"
#include "../../model/model.h"
#include "../dsl/deck_parser.h"
#include "../dsl/file.h"
#include "../writer/writers.h"
#include "commands/register_functions.h"

#include <iostream>
#include <memory>
#include <utility>

namespace fem::io::reader {

Parser::Parser()
    : model_(std::make_shared<model::Model>()), writer_("") {
    configure_documentation_registry();
}

Parser::~Parser() = default;

void Parser::run(const std::string& input_path,
                 const std::string& output_path,
                 const io::writer::WriterFileFormats& writer_formats) {
    model_ = std::make_shared<model::Model>();
    active_loadcase_.reset();
    next_loadcase_id_ = 1;

    io::dsl::Registry registry;
    register_commands(registry);

    io::dsl::File       file(input_path);
    io::dsl::DeckParser deck_parser(registry);
    const io::dsl::Deck deck = deck_parser.parse(file);
    process_deck(deck, input_path, output_path, writer_formats);
}

void Parser::process_deck(const io::dsl::Deck&                 deck,
                          const std::string&                   input_path,
                          const std::string&                   output_path,
                          const io::writer::WriterFileFormats& writer_formats) {
    const auto& root = deck.root();

    root.execute_children("HEADING");

    for (const auto* material : root.children("MATERIAL")) {
        material->enter();
        material->execute_children("ELASTIC");
        material->execute_children("HYPERELASTIC");
        material->execute_children("DENSITY");
        material->execute_children("CONDUCTIVITY");
        material->execute_children("THERMALEXPANSION");
        material->leave();
    }

    root.execute_children("PROFILE");
    root.execute_children("ORIENTATION");
    root.execute_children("AMPLITUDE");

    // Default Part.
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

    // Explicit Parts.
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

    for (const auto* assembly : root.children("ASSEMBLY")) {
        assembly->execute_children("NODE");
        assembly->execute_children("ELEMENT");
    }

    root.execute_children("INSTANCE");
    for (const auto* assembly : root.children("ASSEMBLY")) {
        assembly->execute_children("INSTANCE");
    }

    model_->compile();

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

    model_->build_shell_element_normals();
    root.execute_children("POINTMASS");

    // Root-level Dirichlet and Neumann conditions.
    root.execute_children("SUPPORT");
    root.execute_children("TEMPERATURE");
    root.execute_children("CLOAD");
    root.execute_children("DLOAD");
    root.execute_children("PLOAD");
    root.execute_children("HEATFLUX");
    root.execute_children("CONVECTION");
    root.execute_children("TLOAD");
    root.execute_children("VLOAD");
    root.execute_children("INERTIALLOAD");

    root.execute_children("RBM");
    root.execute_children("CONNECTOR");
    root.execute_children("TIE");
    root.execute_children("CONTACT");
    root.execute_children("EQUATION");

    // Assembly-level Dirichlet and Neumann conditions.
    for (const auto* assembly : root.children("ASSEMBLY")) {
        assembly->execute_children("SUPPORT");
        assembly->execute_children("TEMPERATURE");
        assembly->execute_children("CLOAD");
        assembly->execute_children("DLOAD");
        assembly->execute_children("PLOAD");
        assembly->execute_children("HEATFLUX");
        assembly->execute_children("CONVECTION");
        assembly->execute_children("TLOAD");
        assembly->execute_children("VLOAD");
        assembly->execute_children("INERTIALLOAD");
        assembly->execute_children("RBM");
        assembly->execute_children("CONNECTOR");
        assembly->execute_children("TIE");
        assembly->execute_children("CONTACT");
        assembly->execute_children("EQUATION");
    }

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

    initialize_writers(input_path, output_path, writer_formats);
    for (const auto* loadcase : root.children("LOADCASE")) {
        loadcase->enter();
        loadcase->execute_children();
        loadcase->leave();
    }
    close_writers();
}

void Parser::initialize_writers(const std::string&                   input_path,
                                const std::string&                   output_path,
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

void Parser::close_writers() {
    writer_.close();
}

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

model::Model& Parser::model() {
    logging::error(model_ != nullptr, "Parser: model is not initialized");
    return *model_;
}

const model::Model& Parser::model() const {
    logging::error(model_ != nullptr, "Parser: model is not initialized");
    return *model_;
}

const io::dsl::Registry& Parser::registry() const {
    return documentation_registry_;
}

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

void Parser::end_loadcase() {
    logging::error(active_loadcase_ != nullptr,
        "Parser: cannot end a load case when none is active");

    auto loadcase = std::move(active_loadcase_);
    loadcase->run();
}

loadcase::LoadCase* Parser::active_loadcase() {
    return active_loadcase_.get();
}

void Parser::configure_documentation_registry() {
    documentation_registry_ = io::dsl::Registry{};
    register_commands(documentation_registry_);
}

void Parser::register_commands(io::dsl::Registry& registry) {
    logging::error(model_ != nullptr,
        "Parser: model must exist before registering commands");

    auto& mdl = *model_;

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

    commands::register_heading(registry);
    commands::register_field(registry, mdl);
    commands::register_normal(registry, mdl);
    commands::register_material(registry, mdl);
    commands::register_elastic(registry, mdl);
    commands::register_hyperelastic(registry, mdl);
    commands::register_density(registry, mdl);
    commands::register_conductivity(registry, mdl);
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

    commands::register_cload(registry, mdl);
    commands::register_dload(registry, mdl);
    commands::register_pload(registry, mdl);
    commands::register_heat_flux(registry, mdl);
    commands::register_convection(registry, mdl);
    commands::register_tload(registry, mdl);
    commands::register_vload(registry, mdl);
    commands::register_inertialload(registry, mdl);
    commands::register_rbm(registry, mdl);
    commands::register_support(registry, mdl);
    commands::register_temperature(registry, mdl);
    commands::register_amplitude(registry, mdl);
    commands::register_connector(registry, mdl);
    commands::register_coupling(registry, mdl);
    commands::register_tie(registry, mdl);
    commands::register_contact(registry, mdl);
    commands::register_point_mass(registry, mdl);
    commands::register_overview(registry, mdl);
    commands::register_equation(registry, mdl);

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
