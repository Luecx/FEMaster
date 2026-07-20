#include "parser.h"

#include "../dsl/engine.h"
#include "../dsl/file.h"
#include "../../loadcase/linear_buckling.h"
#include "../../loadcase/linear_eigenfreq.h"
#include "../../loadcase/linear_static.h"
#include "../../loadcase/linear_static_topo.h"
#include "../../loadcase/loadcase.h"
#include "../../model/model.h"
#include "../writer/writers.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <utility>

// Command registration helpers
#include "commands/register_node_count.inl"
#include "commands/register_element_count.inl"
#include "commands/register_surface_count.inl"
#include "commands/register_heading.inl"
#include "commands/register_field.inl"
#include "commands/register_node.inl"
#include "commands/register_nset.inl"
#include "commands/register_elset.inl"
#include "commands/register_surface.inl"
#include "commands/register_sfset.inl"
#include "commands/register_material.inl"
#include "commands/register_elastic.inl"
#include "commands/register_hyperelastic.inl"
#include "commands/register_density.inl"
#include "commands/register_thermalexpansion.inl"
#include "commands/register_cload.inl"
#include "commands/register_dload.inl"
#include "commands/register_pload.inl"
#include "commands/register_tload.inl"
#include "commands/register_vload.inl"
#include "commands/register_inertialload.inl"
#include "commands/register_rbm.inl"
#include "commands/register_support.inl"
#include "commands/register_amplitude.inl"
#include "commands/register_orientation.inl"
#include "commands/register_connector.inl"
#include "commands/register_coupling.inl"
#include "commands/register_tie.inl"
#include "commands/register_contact.inl"
#include "commands/register_profile.inl"
#include "commands/register_solid_section.inl"
#include "commands/register_beam_section.inl"
#include "commands/register_truss_section.inl"
#include "commands/register_shell_section.inl"
#include "commands/register_point_mass.inl"
#include "commands/register_element.inl"
#include "commands/register_overview.inl"
#include "commands/register_loadcase_begin.inl"
#include "commands/register_loadcase_supports.inl"
#include "commands/register_loadcase_loads.inl"
#include "commands/register_loadcase_solver.inl"
#include "commands/register_loadcase_constraintmethod.inl"
#include "commands/register_loadcase_request_stiffness.inl"
#include "commands/register_loadcase_request_stgeom.inl"
#include "commands/register_loadcase_numeigenvalues.inl"
#include "commands/register_loadcase_sigma.inl"
#include "commands/register_loadcase_topodensity.inl"
#include "commands/register_loadcase_topoorient.inl"
#include "commands/register_loadcase_topoexponent.inl"
#include "commands/register_loadcase_constraintsummary.inl"
#include "commands/register_loadcase_nonlinear.inl"

// NEW transient-related commands
#include "commands/register_loadcase_time.inl"
#include "commands/register_loadcase_write_every.inl"
#include "commands/register_loadcase_damping.inl"
#include "commands/register_loadcase_newmark.inl"
#include "commands/register_loadcase_initialvelocity.inl"
#include "commands/register_loadcase_inertiarelief.inl"
#include "commands/register_loadcase_rebalance.inl"

namespace fem::io::reader {

Parser::Parser()
    : m_model(std::make_shared<model::Model>(1, 1, 1)) // tiny placeholder
    , m_writer("")                                      // opened later in run()
{
    // Commands available right away (for documentation mode).
    register_documentation_commands();
}

Parser::~Parser() = default;

// ----------------- Public API -----------------

void Parser::run(const std::string&                  input_path,
                 const std::string&                  output_path,
                 const io::writer::WriterFileFormats& writer_formats) {
    // 1) Count ids needed for ModelData allocation.
    CountData count = run_count_stage(input_path);

    // 2) Rebuild model with correct capacities, reset state
    allocate_model(count);

    // 3) Parse topology and sections before element offsets are built
    run_topology_stage(input_path);
    m_model->assign_sections();
    m_model->_data->initialize_element_enumeration();
    m_model->build_shell_element_normals();

    // 4) Parse all non-topology input and field data.
    run_data_stage(input_path, output_path, writer_formats);
}

void Parser::document(const DocOptions& opts) const {
    using A = DocOptions::Action;
    using F = DocOptions::Format;
    using V = DocOptions::Verbosity;

    // For now, only "text" is implemented. Others can route to same text or be extended later.
    const bool as_text = (opts.format == F::Text);

    if (!as_text) {
        std::cout << "(Note) Only text output is implemented currently. Falling back to text.\n\n";
    }

    switch (opts.action) {
        case A::List:
            m_registry.print_index();
            break;

        case A::Show:
            if (opts.verbosity == V::Compact) m_registry.print_help(opts.cmd, /*compact=*/true);
            else                              m_registry.print_help(opts.cmd, /*compact=*/false);
            break;

        case A::Tokens:
            m_registry.print_tokens(opts.cmd);
            break;

        case A::Variants:
            m_registry.print_variants(opts.cmd);
            break;

        case A::Search:
            m_registry.print_search(opts.query, opts.regex);
            break;

        case A::WhereToken:
            m_registry.print_where_token(opts.query);
            break;

        case A::All:
            m_registry.print_help(/*filter=*/{}, /*compact=*/false);
            break;
    }
}

// ----------------- Accessors -----------------

model::Model& Parser::model() {
    if (!m_model) throw std::runtime_error("Model not initialized.");
    return *m_model;
}
const model::Model& Parser::model() const {
    if (!m_model) throw std::runtime_error("Model not initialized.");
    return *m_model;
}
io::writer::ResultWriters& Parser::writer() { return m_writer; }
const io::writer::ResultWriters& Parser::writer() const { return m_writer; }
io::dsl::Registry& Parser::registry() { return m_registry; }
const io::dsl::Registry& Parser::registry() const { return m_registry; }

// ----------------- Loadcase bookkeeping -----------------

int Parser::next_loadcase_id() { return m_next_loadcase_id++; }

void Parser::set_active_loadcase(std::unique_ptr<loadcase::LoadCase> lc, std::string type) {
    m_active_loadcase      = std::move(lc);
    m_active_loadcase_type = std::move(type);
}
loadcase::LoadCase* Parser::active_loadcase() { return m_active_loadcase.get(); }
const loadcase::LoadCase* Parser::active_loadcase() const { return m_active_loadcase.get(); }
void Parser::clear_active_loadcase() { m_active_loadcase.reset(); m_active_loadcase_type.clear(); }
const std::string& Parser::active_loadcase_type() const {
    static const std::string empty;
    return m_active_loadcase_type.empty() ? empty : m_active_loadcase_type;
}

// ----------------- Parser stages -----------------

Parser::CountData Parser::run_count_stage(const std::string& input_path) {
    CountData count;
    io::dsl::Registry registry;
    register_count_commands(registry, count);
    register_set_commands(registry);
    register_analysis_commands(registry);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("NODE", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELEMENT", io::dsl::ActiveMode::Active);
    registry.set_active_mode("SURFACE", io::dsl::ActiveMode::Active);

    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);

    return count;
}

void Parser::run_topology_stage(const std::string& input_path) {
    io::dsl::Registry registry;
    register_topology_commands(registry);
    register_analysis_commands(registry);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("NODE", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELEMENT", io::dsl::ActiveMode::Active);
    registry.set_active_mode("NSET", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELSET", io::dsl::ActiveMode::Active);
    registry.set_active_mode("SURFACE", io::dsl::ActiveMode::Active);
    registry.set_active_mode("SFSET", io::dsl::ActiveMode::Active);
    registry.set_active_mode("MATERIAL", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELASTIC", io::dsl::ActiveMode::Active);
    registry.set_active_mode("DENSITY", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ORIENTATION", io::dsl::ActiveMode::Active);
    registry.set_active_mode("SHELLSECTION", io::dsl::ActiveMode::Active);

    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);
}

// ----------------- Model & registration -----------------

void Parser::allocate_model(const CountData& count) {
    const int n_nodes    = count.highest_node_id    + 1;
    const int n_elems    = count.highest_element_id + 1;
    const int n_surfaces = count.highest_surface_id + 1;

    m_model = std::make_shared<model::Model>(n_nodes, n_elems, n_surfaces);

    // Reset state bound to an old model
    m_active_loadcase.reset();
    m_active_loadcase_type.clear();
    m_next_loadcase_id = 1;
}

void Parser::run_data_stage(const std::string&                  input_path,
                            const std::string&                  output_path,
                            const io::writer::WriterFileFormats& writer_formats) {
    std::string writer_base = output_path.empty() ? input_path : output_path;
    for (const std::string& ext : {std::string(".res"), std::string(".frd"), std::string(".inp")}) {
        if (writer_base.size() >= ext.size() &&
            writer_base.compare(writer_base.size() - ext.size(), ext.size(), ext) == 0) {
            writer_base.resize(writer_base.size() - ext.size());
            break;
        }
    }

    m_writer = io::writer::ResultWriters(writer_base, writer_formats);
    m_writer.write_model_data(*m_model->_data);

    io::dsl::Registry registry;
    register_topology_commands(registry);
    register_analysis_commands(registry);

    registry.set_active_mode(io::dsl::ActiveMode::Active);
    registry.set_active_mode("NODE", io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("ELEMENT", io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("NSET", io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("ELSET", io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("SURFACE", io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("SFSET", io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("MATERIAL", io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("ELASTIC", io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("DENSITY", io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("ORIENTATION", io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("SHELLSECTION", io::dsl::ActiveMode::ConsumeOnly);

    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);

    m_writer.close();
}

void Parser::register_documentation_commands() {
    m_registry = io::dsl::Registry{};
    register_topology_commands(m_registry);
    register_analysis_commands(m_registry);
    m_registry.set_active_mode(io::dsl::ActiveMode::Active);
}

void Parser::register_count_commands(io::dsl::Registry& reg, CountData& count) {
    commands::register_node_count(reg, [&count](ID id) {
        count.highest_node_id = std::max(count.highest_node_id, static_cast<int>(id));
    });
    commands::register_element_count(reg, [&count](ID id) {
        count.highest_element_id = std::max(count.highest_element_id, static_cast<int>(id));
    });
    commands::register_surface_count(reg, [&count](ID id) {
        count.highest_surface_id = std::max(count.highest_surface_id, static_cast<int>(id));
    });
}

void Parser::register_set_commands(io::dsl::Registry& reg) {
    if (!m_model) throw std::runtime_error("Model must exist before registering commands");

    auto& mdl = *m_model;
    commands::register_nset(reg, mdl);
    commands::register_elset(reg, mdl);
    commands::register_sfset(reg, mdl);
}

void Parser::register_topology_commands(io::dsl::Registry& reg) {
    if (!m_model) throw std::runtime_error("Model must exist before registering commands");

    auto& mdl = *m_model;
    commands::register_node(reg, mdl);
    commands::register_element(reg, mdl);
    commands::register_nset(reg, mdl);
    commands::register_elset(reg, mdl);
    commands::register_surface(reg, mdl);
    commands::register_sfset(reg, mdl);
}

void Parser::register_analysis_commands(io::dsl::Registry& reg) {
    if (!m_model) throw std::runtime_error("Model must exist before registering commands");

    auto& mdl = *m_model;

    // Structural commands
    reg.command("MODEL", [](io::dsl::Command& command) {
        command.allow_if(io::dsl::Condition::parent_is("ROOT"));
        command.doc("Begin a model definition.");
        command.keyword(io::dsl::KeywordSpec::make().key("NAME").optional());
        command.variant(io::dsl::Variant::make());
    });
    reg.command("END", [](io::dsl::Command& command) {
        command.allow_if(io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("End the current load case definition.");
        command.variant(io::dsl::Variant::make());
    });

    // Base model/sections/materials
    commands::register_heading(reg);
    commands::register_field(reg, mdl);
    commands::register_material(reg, mdl);
    commands::register_elastic(reg, mdl);
    commands::register_hyperelastic(reg, mdl);
    commands::register_density(reg, mdl);
    commands::register_thermal_expansion(reg, mdl);

    // Loads & BCs
    commands::register_cload(reg, mdl);
    commands::register_dload(reg, mdl);
    commands::register_pload(reg, mdl);
    commands::register_tload(reg, mdl);
    commands::register_vload(reg, mdl);
    commands::register_inertialload(reg, mdl);
    commands::register_rbm(reg, mdl);
    commands::register_support(reg, mdl);
    commands::register_amplitude(reg, mdl);

    // Orientations & connectors/constraints
    commands::register_orientation(reg, mdl);
    commands::register_connector(reg, mdl);
    commands::register_coupling(reg, mdl);
    commands::register_tie(reg, mdl);
    commands::register_contact(reg, mdl);

    // Profiles & sections & elements
    commands::register_profile(reg, mdl);
    commands::register_solid_section(reg, mdl);
    commands::register_beam_section(reg, mdl);
    commands::register_truss_section(reg, mdl);
    commands::register_shell_section(reg, mdl);
    commands::register_point_mass(reg, mdl);
    commands::register_overview(reg, mdl);

    // Loadcase scaffold
    commands::register_loadcase_begin(reg, *this);
    commands::register_loadcase_supports(reg, *this);
    commands::register_loadcase_loads(reg, *this);
    commands::register_loadcase_solver(reg, *this);
    commands::register_loadcase_constraintmethod(reg, *this);
    commands::register_loadcase_request_stiffness(reg, *this);
    commands::register_loadcase_request_stgeom(reg, *this);
    commands::register_loadcase_numeigenvalues(reg, *this);
    commands::register_loadcase_sigma(reg, *this);
    commands::register_loadcase_topodensity(reg, *this);
    commands::register_loadcase_topoorient(reg, *this);
    commands::register_loadcase_topoexponent(reg, *this);
    commands::register_loadcase_constraintsummary(reg, *this);
    commands::register_loadcase_nonlinear(reg, *this);

    // NEW: Transient-specific loadcase commands
    commands::register_loadcase_time(reg, *this);
    commands::register_loadcase_write_every(reg, *this);
    commands::register_loadcase_damping(reg, *this);
    commands::register_loadcase_newmark(reg, *this);
    commands::register_loadcase_initialvelocity(reg, *this);
    commands::register_loadcase_inertiarelief(reg, *this);
    commands::register_loadcase_rebalance(reg, *this);
}
} // namespace fem::io::reader
