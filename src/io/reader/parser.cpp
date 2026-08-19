/**
 * @file parser.cpp
 * @brief Implements dependency-ordered FEMaster input-deck parsing.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "parser.h"

#include "../dsl/engine.h"
#include "../dsl/file.h"
#include "../writer/writers.h"
#include "../../core/logging.h"
#include "../../loadcase/loadcase.h"
#include "../../model/model.h"

#include <iostream>
#include <memory>
#include <utility>

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
#include "commands/register_field.inl"
#include "commands/register_heading.inl"
#include "commands/register_hyperelastic.inl"
#include "commands/register_inertialload.inl"
#include "commands/register_instance.inl"
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
#include "commands/register_sfset.inl"
#include "commands/register_shell_section.inl"
#include "commands/register_solid_section.inl"
#include "commands/register_support.inl"
#include "commands/register_surface.inl"
#include "commands/register_thermalexpansion.inl"
#include "commands/register_tie.inl"
#include "commands/register_tload.inl"
#include "commands/register_truss_section.inl"
#include "commands/register_vload.inl"

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

namespace fem::io::reader {

Parser::Parser()
    : m_model(std::make_shared<model::Model>()),
      m_writer("") {
    register_documentation_commands();
}

Parser::~Parser() = default;

void Parser::run(const std::string& input_path,
                 const std::string& output_path,
                 const io::writer::WriterFileFormats& writer_formats) {
    m_model = std::make_shared<model::Model>();
    m_active_loadcase.reset();
    m_active_loadcase_type.clear();
    m_next_loadcase_id = 1;

    run_definition_stage(input_path);
    run_topology_stage(input_path);
    m_model->compile();
    run_field_stage(input_path);
    m_model->build_shell_element_normals();
    run_data_stage(input_path, output_path, writer_formats);
}

void Parser::document(const DocOptions& opts) const {
    using A = DocOptions::Action;
    using F = DocOptions::Format;
    using V = DocOptions::Verbosity;

    if (opts.format != F::Text) {
        std::cout << "(Note) Only text output is implemented currently. Falling back to text.\n\n";
    }

    switch (opts.action) {
        case A::List:       m_registry.print_index(); break;
        case A::Show:       m_registry.print_help(opts.cmd, opts.verbosity == V::Compact); break;
        case A::Tokens:     m_registry.print_tokens(opts.cmd); break;
        case A::Variants:   m_registry.print_variants(opts.cmd); break;
        case A::Search:     m_registry.print_search(opts.query, opts.regex); break;
        case A::WhereToken: m_registry.print_where_token(opts.query); break;
        case A::All:        m_registry.print_help({}, false); break;
    }
}

model::Model& Parser::model() {
    logging::error(m_model != nullptr,
        "Parser: model is not initialized");
    return *m_model;
}

const model::Model& Parser::model() const {
    logging::error(m_model != nullptr,
        "Parser: model is not initialized");
    return *m_model;
}

io::writer::ResultWriters& Parser::writer() { return m_writer; }
const io::writer::ResultWriters& Parser::writer() const { return m_writer; }
io::dsl::Registry& Parser::registry() { return m_registry; }
const io::dsl::Registry& Parser::registry() const { return m_registry; }

int Parser::next_loadcase_id() { return m_next_loadcase_id++; }

void Parser::set_active_loadcase(std::unique_ptr<loadcase::LoadCase> loadcase, std::string type) {
    m_active_loadcase      = std::move(loadcase);
    m_active_loadcase_type = std::move(type);
}

loadcase::LoadCase* Parser::active_loadcase() { return m_active_loadcase.get(); }
const loadcase::LoadCase* Parser::active_loadcase() const { return m_active_loadcase.get(); }

void Parser::clear_active_loadcase() {
    m_active_loadcase.reset();
    m_active_loadcase_type.clear();
}

const std::string& Parser::active_loadcase_type() const {
    static const std::string empty;
    return m_active_loadcase_type.empty() ? empty : m_active_loadcase_type;
}

void Parser::run_definition_stage(const std::string& input_path) {
    io::dsl::Registry registry;
    configure_definition_stage(registry);
    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);
}

void Parser::run_topology_stage(const std::string& input_path) {
    io::dsl::Registry registry;
    configure_topology_stage(registry);
    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);
}

void Parser::run_field_stage(const std::string& input_path) {
    io::dsl::Registry registry;
    configure_field_stage(registry);
    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);
}

void Parser::run_data_stage(const std::string& input_path,
                            const std::string& output_path,
                            const io::writer::WriterFileFormats& writer_formats) {
    std::string writer_base = output_path.empty() ? input_path : output_path;
    for (const std::string& ext : {std::string(".res"), std::string(".frd"), std::string(".inp")}) {
        if (writer_base.size() >= ext.size()
         && writer_base.compare(writer_base.size() - ext.size(), ext.size(), ext) == 0) {
            writer_base.resize(writer_base.size() - ext.size());
            break;
        }
    }

    m_writer = io::writer::ResultWriters(writer_base, writer_formats);
    m_writer.write_model_data(*m_model->_data);

    io::dsl::Registry registry;
    configure_data_stage(registry);
    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);
    m_writer.close();
}

void Parser::configure_definition_stage(io::dsl::Registry& registry) {
    register_topology_commands(registry);
    register_analysis_commands(registry);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("MATERIAL", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELASTIC", io::dsl::ActiveMode::Active);
    registry.set_active_mode("HYPERELASTIC", io::dsl::ActiveMode::Active);
    registry.set_active_mode("DENSITY", io::dsl::ActiveMode::Active);
    registry.set_active_mode("THERMALEXPANSION", io::dsl::ActiveMode::Active);
    registry.set_active_mode("PROFILE", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ORIENTATION", io::dsl::ActiveMode::Active);
}

void Parser::configure_topology_stage(io::dsl::Registry& registry) {
    register_topology_commands(registry);
    register_analysis_commands(registry);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    for (const std::string& command : {
        "PART", "ENDPART", "ASSEMBLY", "ENDASSEMBLY", "INSTANCE", "ENDINSTANCE",
        "NODE", "ELEMENT", "NSET", "ELSET", "SURFACE", "SFSET",
        "SOLIDSECTION", "BEAMSECTION", "TRUSSSECTION", "SHELLSECTION"
    }) {
        registry.set_active_mode(command, io::dsl::ActiveMode::Active);
    }
}

void Parser::configure_field_stage(io::dsl::Registry& registry) {
    register_topology_commands(registry);
    register_analysis_commands(registry);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    for (const std::string& command : {
        "ASSEMBLY", "ENDASSEMBLY", "NSET", "ELSET", "SURFACE", "SFSET", "FIELD", "NORMAL"
    }) {
        registry.set_active_mode(command, io::dsl::ActiveMode::Active);
    }
}

void Parser::configure_data_stage(io::dsl::Registry& registry) {
    register_topology_commands(registry);
    register_analysis_commands(registry);
    registry.set_active_mode(io::dsl::ActiveMode::Active);

    for (const std::string& command : {
        "PART", "ENDPART", "ASSEMBLY", "ENDASSEMBLY", "INSTANCE", "ENDINSTANCE",
        "NODE", "ELEMENT", "NSET", "ELSET", "SURFACE", "SFSET",
        "MATERIAL", "ELASTIC", "HYPERELASTIC", "DENSITY", "THERMALEXPANSION",
        "PROFILE", "ORIENTATION", "SOLIDSECTION", "BEAMSECTION", "TRUSSSECTION",
        "SHELLSECTION", "FIELD", "NORMAL"
    }) {
        registry.set_active_mode(command, io::dsl::ActiveMode::ConsumeOnly);
    }
}

void Parser::register_documentation_commands() {
    m_registry = io::dsl::Registry{};
    register_topology_commands(m_registry);
    register_analysis_commands(m_registry);
    m_registry.set_active_mode(io::dsl::ActiveMode::Active);
}

void Parser::register_topology_commands(io::dsl::Registry& registry) {
    logging::error(m_model != nullptr,
        "Parser: model must exist before registering topology commands");

    auto& mdl = *m_model;
    auto assembly_scope = std::make_shared<bool>(false);

    commands::register_part        (registry, mdl);
    commands::register_end_part    (registry, mdl);
    commands::register_assembly    (registry, assembly_scope);
    commands::register_end_assembly(registry, assembly_scope);
    commands::register_instance    (registry, mdl);
    commands::register_end_instance(registry);
    commands::register_node        (registry, mdl, assembly_scope);
    commands::register_element     (registry, mdl, assembly_scope);
    commands::register_nset        (registry, mdl, assembly_scope);
    commands::register_elset       (registry, mdl, assembly_scope);
    commands::register_surface     (registry, mdl, assembly_scope);
    commands::register_sfset       (registry, mdl, assembly_scope);
}

void Parser::register_analysis_commands(io::dsl::Registry& registry) {
    logging::error(m_model != nullptr,
        "Parser: model must exist before registering analysis commands");

    auto& mdl = *m_model;

    registry.command("MODEL", [](io::dsl::Command& command) {
        command.allow_if(io::dsl::Condition::parent_is("ROOT"));
        command.keyword(io::dsl::KeywordSpec::make().key("NAME").optional());
        command.variant(io::dsl::Variant::make());
    });
    registry.command("END", [](io::dsl::Command& command) {
        command.allow_if(io::dsl::Condition::parent_is("LOADCASE"));
        command.variant(io::dsl::Variant::make());
    });

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
