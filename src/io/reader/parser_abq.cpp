/**
 * @file parser_abq.cpp
 * @brief Implements staged Abaqus input-deck parsing.
 *
 * The initial reader supports only HEADING, NODE, ELEMENT, NSET and ELSET.
 * HEADING, NODE, NSET and ELSET reuse the native FEMaster command definitions;
 * ELEMENT has a dedicated Abaqus registration for Abaqus element type names.
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#include "parser_abq.h"

#include "../dsl/engine.h"
#include "../dsl/file.h"
#include "../dsl/registry.h"
#include "../../model/model.h"
#include "../writer/writers.h"

#include <algorithm>
#include <string>

#include "commands/register_node_count.inl"
#include "commands/register_heading.inl"
#include "commands/register_node.inl"
#include "commands/register_nset.inl"
#include "commands/register_elset.inl"
#include "commands_abq/register_element.inl"

namespace fem::io::reader {

ParserAbq::ParserAbq()
    : m_model(std::make_shared<model::Model>(1, 1, 1))
    , m_writer("") {}

ParserAbq::~ParserAbq() = default;

void ParserAbq::run(const std::string&                   input_path,
                    const std::string&                   output_path,
                    const io::writer::WriterFileFormats& writer_formats) {
    CountData count = run_count_stage(input_path);
    allocate_model(count);

    run_topology_stage(input_path);
    m_model->assign_sections();
    m_model->_data->initialize_element_enumeration();

    run_field_stage(input_path);
    m_model->build_shell_element_normals();

    run_data_stage(input_path, output_path, writer_formats);
}

ParserAbq::CountData ParserAbq::run_count_stage(const std::string& input_path) {
    CountData count;
    io::dsl::Registry registry;

    commands::register_heading(registry);
    commands::register_node_count(registry, [&count](ID id) {
        count.highest_node_id = std::max(count.highest_node_id, static_cast<int>(id));
    });
    commands_abq::register_element_count(registry, [&count](ID id) {
        count.highest_element_id = std::max(count.highest_element_id, static_cast<int>(id));
    });
    commands::register_nset(registry, *m_model);
    commands::register_elset(registry, *m_model);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("NODE", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELEMENT", io::dsl::ActiveMode::Active);

    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);

    return count;
}

void ParserAbq::run_topology_stage(const std::string& input_path) {
    io::dsl::Registry registry;

    commands::register_heading(registry);
    commands::register_node(registry, *m_model);
    commands_abq::register_element(registry, *m_model);
    commands::register_nset(registry, *m_model);
    commands::register_elset(registry, *m_model);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);
    registry.set_active_mode("NODE", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELEMENT", io::dsl::ActiveMode::Active);
    registry.set_active_mode("NSET", io::dsl::ActiveMode::Active);
    registry.set_active_mode("ELSET", io::dsl::ActiveMode::Active);

    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);
}

void ParserAbq::run_field_stage(const std::string& input_path) {
    io::dsl::Registry registry;

    commands::register_heading(registry);
    commands::register_node(registry, *m_model);
    commands_abq::register_element(registry, *m_model);
    commands::register_nset(registry, *m_model);
    commands::register_elset(registry, *m_model);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);

    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);
}

void ParserAbq::run_data_stage(const std::string&                   input_path,
                               const std::string&                   output_path,
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
    commands::register_heading(registry);
    commands::register_node(registry, *m_model);
    commands_abq::register_element(registry, *m_model);
    commands::register_nset(registry, *m_model);
    commands::register_elset(registry, *m_model);

    registry.set_active_mode(io::dsl::ActiveMode::ConsumeOnly);

    io::dsl::File file(input_path);
    io::dsl::Engine engine(registry);
    engine.run(file);

    m_writer.close();
}

void ParserAbq::allocate_model(const CountData& count) {
    const int n_nodes    = count.highest_node_id + 1;
    const int n_elems    = count.highest_element_id + 1;
    const int n_surfaces = count.highest_surface_id + 1;
    m_model = std::make_shared<model::Model>(n_nodes, n_elems, n_surfaces);
}

} // namespace fem::io::reader
