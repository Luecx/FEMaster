#include "parser.h"

#include <algorithm>
#include <cstddef>
#include <stdexcept>

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
#include "commands/register_nset.inl"
#include "commands/register_elset.inl"
#include "commands/register_material.inl"
#include "commands/register_elastic.inl"
#include "commands/register_density.inl"
#include "commands/register_thermalexpansion.inl"
#include "commands/register_cload.inl"
#include "commands/register_dload.inl"
#include "commands/register_pload.inl"
#include "commands/register_tload.inl"
#include "commands/register_vload.inl"
#include "commands/register_support.inl"
#include "commands/register_temperature.inl"
#include "commands/register_orientation.inl"
#include "commands/register_connector.inl"
#include "commands/register_coupling.inl"
#include "commands/register_tie.inl"
#include "commands/register_profile.inl"
#include "commands/register_solid_section.inl"
#include "commands/register_beam_section.inl"
#include "commands/register_shell_section.inl"
#include "commands/register_point_mass_section.inl"
#include "commands/register_element.inl"
#include "commands/register_loadcase_begin.inl"
#include "commands/register_loadcase_supports.inl"
#include "commands/register_loadcase_loads.inl"
#include "commands/register_loadcase_solver.inl"
#include "commands/register_loadcase_request_stiffness.inl"
#include "commands/register_loadcase_request_stgeom.inl"
#include "commands/register_loadcase_numeigenvalues.inl"
#include "commands/register_loadcase_sigma.inl"
#include "commands/register_loadcase_topodensity.inl"
#include "commands/register_loadcase_topoorient.inl"
#include "commands/register_loadcase_topoexponent.inl"
#include "commands/register_loadcase_constraintsummary.inl"

namespace fem::input_decks {

Parser::Parser(std::string input_path, std::string output_path)
    : m_input_path(std::move(input_path)),
      m_output_path(std::move(output_path)),
      m_writer(m_output_path) {}

Parser::~Parser() = default;

void Parser::preprocess() {
    m_analyse = {};

    dsl::File file(m_input_path);
    dsl::Line current;
    current = file.next_line();

    while (current.type() != dsl::LineType::END_OF_FILE) {
        if (current.type() != dsl::LineType::KEYWORD_LINE) {
            current = file.next_line();
            continue;
        }

        const std::string cmd = current.command();
        if (cmd == "NODE") {
            analyse_nodes(file, current);
            continue;
        }
        if (cmd == "ELEMENT") {
            analyse_elements(file, current);
            continue;
        }
        if (cmd == "SURFACE") {
            analyse_surfaces(file, current);
            continue;
        }

        current = file.next_line();
    }

    m_model = std::make_shared<model::Model>(m_analyse.highest_node_id + 1,
                                             m_analyse.highest_element_id + 1,
                                             m_analyse.highest_surface_id + 1);
    m_commands_registered = false;
    m_next_loadcase_id = 1;
    m_active_loadcase.reset();
    m_active_loadcase_type.clear();
}

void Parser::parse() {
    if (!m_model) {
        preprocess();
    }

    m_writer.open(m_output_path);
    register_commands();

    dsl::File file(m_input_path);
    dsl::Engine engine(m_registry);
    engine.run(file);

    m_writer.close();
}

model::Model& Parser::model() {
    if (!m_model) {
        throw std::runtime_error("Model not initialized. Call preprocess() first.");
    }
    return *m_model;
}

const model::Model& Parser::model() const {
    if (!m_model) {
        throw std::runtime_error("Model not initialized. Call preprocess() first.");
    }
    return *m_model;
}

reader::Writer& Parser::writer() {
    return m_writer;
}

const reader::Writer& Parser::writer() const {
    return m_writer;
}

dsl::Registry& Parser::registry() {
    return m_registry;
}

void Parser::register_all_commands() {
    if (!m_model) {
        throw std::runtime_error("Model not initialized. Call preprocess() first.");
    }
    register_commands();
}

void Parser::prepare_for_documentation() {
    if (!m_model) {
        m_model = std::make_shared<model::Model>(1, 1, 1);
        m_commands_registered = false;
        m_next_loadcase_id = 1;
        m_active_loadcase.reset();
        m_active_loadcase_type.clear();
    }
}

int Parser::next_loadcase_id() {
    return m_next_loadcase_id++;
}

void Parser::set_active_loadcase(std::unique_ptr<loadcase::LoadCase> lc, std::string type) {
    m_active_loadcase = std::move(lc);
    m_active_loadcase_type = std::move(type);
}

loadcase::LoadCase* Parser::active_loadcase() {
    return m_active_loadcase.get();
}

const loadcase::LoadCase* Parser::active_loadcase() const {
    return m_active_loadcase.get();
}

void Parser::clear_active_loadcase() {
    m_active_loadcase.reset();
    m_active_loadcase_type.clear();
}

const std::string& Parser::active_loadcase_type() const {
    static const std::string empty;
    return m_active_loadcase_type.empty() ? empty : m_active_loadcase_type;
}

void Parser::analyse_nodes(dsl::File& file, dsl::Line& current) {
    while (true) {
        current = file.next_line();
        if (current.type() != dsl::LineType::DATA_LINE) {
            break;
        }
        int node_id = std::stoi(current.values()[0]);
        m_analyse.highest_node_id = std::max(m_analyse.highest_node_id, node_id);
    }
}

static std::size_t required_nodes_for_type(const std::string& type) {
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
    if (type == "S4" || type == "MITC4") return 4;
    if (type == "S6") return 6;
    if (type == "S8") return 8;
    if (type == "P") return 1;
    return 0;
}

void Parser::analyse_elements(dsl::File& file, dsl::Line& current) {
    const dsl::Keys keys = dsl::Keys::from_keyword_line(current);
    const std::string type = keys.raw("TYPE");
    const std::size_t req = required_nodes_for_type(type);

    if (req == 0) {
        throw std::runtime_error("Unknown element type during preprocessing: " + type);
    }

    while (true) {
        current = file.next_line();
        if (current.type() != dsl::LineType::DATA_LINE) {
            break;
        }

        int element_id = std::stoi(current.values()[0]);
        m_analyse.highest_element_id = std::max(m_analyse.highest_element_id, element_id);

        std::size_t have = current.values().size() - 1; // minus element id
        while (have < req) {
            current = file.next_line();
            if (current.type() != dsl::LineType::DATA_LINE) {
                throw std::runtime_error("Unexpected end of element connectivity while preprocessing");
            }
            have += current.values().size();
        }
    }
}

void Parser::analyse_surfaces(dsl::File& file, dsl::Line& current) {
    while (true) {
        current = file.next_line();
        if (current.type() != dsl::LineType::DATA_LINE) {
            break;
        }
        try {
            int surface_id = std::stoi(current.values()[0]);
            m_analyse.highest_surface_id = std::max(m_analyse.highest_surface_id, surface_id);
        } catch (const std::invalid_argument&) {
        }
    }
}

void Parser::register_commands() {
    if (m_commands_registered) {
        return;
    }

    if (!m_model) {
        throw std::runtime_error("Model must exist before registering commands");
    }

    auto& reg = m_registry;
    auto& mdl = *m_model;

    commands::register_node(reg, mdl);
    commands::register_nset(reg, mdl);
    commands::register_elset(reg, mdl);
    commands::register_material(reg, mdl);
    commands::register_elastic(reg, mdl);
    commands::register_density(reg, mdl);
    commands::register_thermal_expansion(reg, mdl);
    commands::register_cload(reg, mdl);
    commands::register_dload(reg, mdl);
    commands::register_pload(reg, mdl);
    commands::register_tload(reg, mdl);
    commands::register_vload(reg, mdl);
    commands::register_support(reg, mdl);
    commands::register_temperature(reg, mdl);
    commands::register_orientation(reg, mdl);
    commands::register_connector(reg, mdl);
    commands::register_coupling(reg, mdl);
    commands::register_tie(reg, mdl);
    commands::register_profile(reg, mdl);
    commands::register_solid_section(reg, mdl);
    commands::register_beam_section(reg, mdl);
    commands::register_shell_section(reg, mdl);
    commands::register_point_mass_section(reg, mdl);
    commands::register_element(reg, mdl);

    commands::register_loadcase_begin(reg, *this);
    commands::register_loadcase_supports(reg, *this);
    commands::register_loadcase_loads(reg, *this);
    commands::register_loadcase_solver(reg, *this);
    commands::register_loadcase_request_stiffness(reg, *this);
    commands::register_loadcase_request_stgeom(reg, *this);
    commands::register_loadcase_numeigenvalues(reg, *this);
    commands::register_loadcase_sigma(reg, *this);
    commands::register_loadcase_topodensity(reg, *this);
    commands::register_loadcase_topoorient(reg, *this);
    commands::register_loadcase_topoexponent(reg, *this);
    commands::register_loadcase_constraintsummary(reg, *this);

    m_commands_registered = true;
}

} // namespace fem::input_decks
