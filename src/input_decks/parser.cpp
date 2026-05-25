#include "parser.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <iostream>
#include <utility>

#include "../core/logging.h"
#include "../dsl/engine.h"
#include "../dsl/file.h"
#include "../dsl/keys.h"
#include "../loadcase/loadcase.h"
#include "../loadcase/linear_buckling.h"
#include "../loadcase/linear_eigenfreq.h"
#include "../loadcase/linear_static.h"
#include "../loadcase/linear_static_topo.h"
#include "../model/model.h"
#include "../reader/writer.h"

// Command registration helpers
#include "commands/register_node.inl"
#include "commands/register_field.inl"
#include "commands/register_nset.inl"
#include "commands/register_elset.inl"
#include "commands/register_surface.inl"
#include "commands/register_sfset.inl"
#include "commands/register_material.inl"
#include "commands/register_elastic.inl"
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

// NEW transient-related commands
#include "commands/register_loadcase_time.inl"
#include "commands/register_loadcase_write_every.inl"
#include "commands/register_loadcase_damping.inl"
#include "commands/register_loadcase_newmark.inl"
#include "commands/register_loadcase_initialvelocity.inl"
#include "commands/register_loadcase_inertiarelief.inl"
#include "commands/register_loadcase_rebalance.inl"

namespace fem::input_decks {

Parser::Parser()
    : m_model(std::make_shared<model::Model>(1, 1, 1)) // tiny placeholder
    , m_writer("")                                      // opened later in run()
{
    // Commands available right away (for documentation mode).
    register_commands();
}

Parser::~Parser() = default;

// ----------------- Public API -----------------

void Parser::run(const std::string& input_path, const std::string& output_path) {
    // 1) Prepass: learn topology and capacities
    AnalyseData A = preprocess_ids(input_path);

    // 2) Rebuild model with correct capacities, reset state
    build_model_from(A);

    // 3) Populate topology and build the element-dependent field enumerations once.
    replay_topology(A);
    m_model->_data->initialize_element_enumeration();

    // 4) Re-register commands so all closures capture the NEW model
    register_commands();

    // 5) Parse all non-topology input and field data
    m_writer.open(output_path);

    dsl::File file(input_path);
    dsl::Engine engine(m_registry);
    engine.run(file);

    m_writer.close();
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
reader::Writer& Parser::writer() { return m_writer; }
const reader::Writer& Parser::writer() const { return m_writer; }
dsl::Registry& Parser::registry() { return m_registry; }
const dsl::Registry& Parser::registry() const { return m_registry; }

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

// ----------------- Preprocess (ID scan) -----------------

Parser::AnalyseData Parser::preprocess_ids(const std::string& input_path) {
    AnalyseData A{};
    dsl::File file(input_path);
    dsl::Line current;
    current = file.next_line();

    while (current.type() != dsl::LineType::END_OF_FILE) {
        if (current.type() != dsl::LineType::KEYWORD_LINE) {
            current = file.next_line();
            continue;
        }

        const std::string cmd = current.command();
        if (cmd == "NODE")     { analyse_nodes(file, current, A);     continue; }
        if (cmd == "ELEMENT")  { analyse_elements(file, current, A);  continue; }
        if (cmd == "SURFACE")  { analyse_surfaces(file, current, A);  continue; }

        current = file.next_line();
    }

    return A;
}

std::size_t Parser::required_nodes_for_type(const std::string& type) {
    if (type == "C3D4") return 4;
    if (type == "C3D5") return 5;
    if (type == "C3D6") return 6;
    if (type == "C3D8") return 8;
    if (type == "C3D10") return 10;
    if (type == "C3D15") return 15;
    if (type == "C3D20" || type == "C3D20R") return 20;
    if (type == "B33") return 2;
    if (type == "T3") return 2;
    if (type == "S3") return 3;
    if (type == "S4" || type == "MITC4" || type == "QSPT") return 4;
    if (type == "S6") return 6;
    if (type == "S8") return 8;
    if (type == "P")  return 1;
    return 0;
}

void Parser::analyse_nodes(dsl::File& file, dsl::Line& current, AnalyseData& A) {
    const dsl::Keys keys = dsl::Keys::from_keyword_line(current);
    const std::string nset = keys.has("NSET") ? keys.raw("NSET") : std::string(SET_NODE_ALL);

    while (true) {
        current = file.next_line();
        if (current.type() != dsl::LineType::DATA_LINE) break;
        const auto& values = current.values();
        int node_id = std::stoi(values[0]);
        A.highest_node_id = std::max(A.highest_node_id, node_id);
        A.nodes.push_back(NodeRecord{
            static_cast<ID>(node_id),
            values.size() > 1 && !values[1].empty() ? static_cast<Precision>(std::stod(values[1])) : Precision{0},
            values.size() > 2 && !values[2].empty() ? static_cast<Precision>(std::stod(values[2])) : Precision{0},
            values.size() > 3 && !values[3].empty() ? static_cast<Precision>(std::stod(values[3])) : Precision{0},
            nset
        });
    }
}

void Parser::analyse_elements(dsl::File& file, dsl::Line& current, AnalyseData& A) {
    const dsl::Keys keys = dsl::Keys::from_keyword_line(current);
    const std::string type = keys.raw("TYPE");
    const std::string elset = keys.has("ELSET") ? keys.raw("ELSET") : std::string(SET_ELEM_ALL);
    const std::size_t req  = required_nodes_for_type(type);

    if (req == 0) {
        throw std::runtime_error("Unknown element type during preprocessing: " + type);
    }

    while (true) {
        current = file.next_line();
        if (current.type() != dsl::LineType::DATA_LINE) break;

        const auto& first_values = current.values();
        int element_id = std::stoi(first_values[0]);
        A.highest_element_id = std::max(A.highest_element_id, element_id);

        ElementRecord record;
        record.id = static_cast<ID>(element_id);
        record.type = type;
        record.elset = elset;
        record.nodes.reserve(req);

        for (std::size_t i = 1; i < first_values.size(); ++i) {
            record.nodes.push_back(static_cast<ID>(std::stoi(first_values[i])));
        }

        std::size_t have = record.nodes.size();
        while (have < req) {
            current = file.next_line();
            if (current.type() != dsl::LineType::DATA_LINE) {
                throw std::runtime_error("Unexpected end of element connectivity while preprocessing");
            }
            for (const auto& value : current.values()) {
                record.nodes.push_back(static_cast<ID>(std::stoi(value)));
            }
            have = record.nodes.size();
        }
        A.elements.push_back(std::move(record));
    }
}

void Parser::analyse_surfaces(dsl::File& file, dsl::Line& current, AnalyseData& A) {
    while (true) {
        current = file.next_line();
        if (current.type() != dsl::LineType::DATA_LINE) break;
        try {
            int surface_id = std::stoi(current.values()[0]);
            A.highest_surface_id = std::max(A.highest_surface_id, surface_id);
        } catch (const std::invalid_argument&) {
        }
    }
}

// ----------------- Model & registration -----------------

void Parser::build_model_from(const AnalyseData& A) {
    const int n_nodes    = A.highest_node_id    + 1;
    const int n_elems    = A.highest_element_id + 1;
    const int n_surfaces = A.highest_surface_id + 1;

    m_model = std::make_shared<model::Model>(n_nodes, n_elems, n_surfaces);

    // Reset state bound to an old model
    m_active_loadcase.reset();
    m_active_loadcase_type.clear();
    m_next_loadcase_id = 1;

    // Fresh registry so all closures capture the new Model references
    m_registry = dsl::Registry{};
    m_commands_registered = false;
}

namespace {
template<class Elem, std::size_t... I>
void replay_regular_element_impl(model::Model& model, ID id, const std::vector<ID>& nodes, std::index_sequence<I...>) {
    model.set_element<Elem>(id, nodes[I]...);
}

template<class Elem, std::size_t N>
void replay_regular_element(model::Model& model, ID id, const std::vector<ID>& nodes) {
    logging::error(nodes.size() == N,
                   "ELEMENT ", id, ": expected ", N, " nodes, got ", nodes.size());
    replay_regular_element_impl<Elem>(model, id, nodes, std::make_index_sequence<N>{});
}
}

void Parser::replay_topology(const AnalyseData& A) {
    auto& mdl = model();

    for (const NodeRecord& node : A.nodes) {
        mdl._data->node_sets.activate(node.nset);
        mdl.set_node(node.id, node.x, node.y, node.z);
    }

    for (const ElementRecord& element : A.elements) {
        mdl._data->elem_sets.activate(element.elset);

        const auto& type = element.type;
        const auto& nodes = element.nodes;
        if (type == "C3D4") replay_regular_element<model::C3D4, 4>(mdl, element.id, nodes);
        else if (type == "C3D6") replay_regular_element<model::C3D6, 6>(mdl, element.id, nodes);
        else if (type == "C3D8") replay_regular_element<model::C3D8, 8>(mdl, element.id, nodes);
        else if (type == "C3D10") replay_regular_element<model::C3D10, 10>(mdl, element.id, nodes);
        else if (type == "C3D15") replay_regular_element<model::C3D15, 15>(mdl, element.id, nodes);
        else if (type == "C3D20") replay_regular_element<model::C3D20, 20>(mdl, element.id, nodes);
        else if (type == "C3D20R") replay_regular_element<model::C3D20R, 20>(mdl, element.id, nodes);
        else if (type == "T3") replay_regular_element<model::T3, 2>(mdl, element.id, nodes);
        else if (type == "S3") replay_regular_element<model::S3, 3>(mdl, element.id, nodes);
        else if (type == "S4") replay_regular_element<model::S4, 4>(mdl, element.id, nodes);
        else if (type == "QSPT") replay_regular_element<model::QSPT, 4>(mdl, element.id, nodes);
        else if (type == "MITC4") replay_regular_element<model::MITC4, 4>(mdl, element.id, nodes);
        else if (type == "S6") replay_regular_element<model::S6, 6>(mdl, element.id, nodes);
        else if (type == "S8") replay_regular_element<model::S8, 8>(mdl, element.id, nodes);
        else if (type == "B33") {
            logging::error(nodes.size() == 2 || nodes.size() == 3,
                           "ELEMENT ", element.id, ": B33 expects 2 or 3 node ids, got ", nodes.size());
            const ID orientation = nodes.size() == 3 ? nodes[2] : ID{-1};
            mdl.set_beam_element<model::B33>(element.id, orientation, nodes[0], nodes[1]);
        } else if (type == "C3D5") {
            logging::error(nodes.size() == 5,
                           "ELEMENT ", element.id, ": C3D5 expects 5 node ids, got ", nodes.size());
            mdl.set_element<model::C3D8>(element.id,
                                         nodes[0], nodes[1], nodes[2], nodes[3],
                                         nodes[4], nodes[4], nodes[4], nodes[4]);
        } else {
            throw std::runtime_error("Unknown element type during topology replay: " + type);
        }
    }
}

void Parser::register_commands() {
    if (m_commands_registered) return;
    if (!m_model) throw std::runtime_error("Model must exist before registering commands");

    auto& reg = m_registry;
    auto& mdl = *m_model;

    // Base model/sections/materials
    commands::register_node(reg, mdl);
    commands::register_field(reg, mdl);
    commands::register_nset(reg, mdl);
    commands::register_elset(reg, mdl);
    commands::register_surface(reg, mdl);
    commands::register_sfset(reg, mdl);
    commands::register_material(reg, mdl);
    commands::register_elastic(reg, mdl);
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

    // Profiles & sections & elements
    commands::register_profile(reg, mdl);
    commands::register_solid_section(reg, mdl);
    commands::register_beam_section(reg, mdl);
    commands::register_truss_section(reg, mdl);
    commands::register_shell_section(reg, mdl);
    commands::register_point_mass(reg, mdl);
    commands::register_element(reg, mdl);
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

    // NEW: Transient-specific loadcase commands
    commands::register_loadcase_time(reg, *this);
    commands::register_loadcase_write_every(reg, *this);
    commands::register_loadcase_damping(reg, *this);
    commands::register_loadcase_newmark(reg, *this);
    commands::register_loadcase_initialvelocity(reg, *this);
    commands::register_loadcase_inertiarelief(reg, *this);
    commands::register_loadcase_rebalance(reg, *this);

    m_commands_registered = true;
}
} // namespace fem::input_decks
