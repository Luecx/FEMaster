/**
 * @file parser_abq.h
 * @brief Declares the staged Abaqus input-deck parser.
 *
 * The Abaqus reader follows the same four-pass model-construction order as the
 * native FEMaster reader. Only a deliberately small initial keyword subset is
 * registered; unsupported Abaqus keywords remain hard errors.
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include "../writer/writers.h"

#include <memory>
#include <string>

namespace fem {
namespace model { struct Model; }

namespace io::reader {

/**
 * @class ParserAbq
 * @brief Executes the staged Abaqus input-deck workflow.
 */
class ParserAbq {
public:
    ParserAbq();
    ~ParserAbq();

    void run(const std::string&                   input_path,
             const std::string&                   output_path,
             const io::writer::WriterFileFormats& writer_formats = io::writer::WriterFileFormats());

private:
    struct CountData {
        int highest_node_id    = -1;
        int highest_element_id = -1;
        int highest_surface_id = -1;
    };

    CountData run_count_stage(const std::string& input_path);
    void run_topology_stage(const std::string& input_path);
    void run_field_stage(const std::string& input_path);
    void run_data_stage(const std::string&                   input_path,
                        const std::string&                   output_path,
                        const io::writer::WriterFileFormats& writer_formats);

    void allocate_model(const CountData& count);

private:
    std::shared_ptr<model::Model> m_model;
    io::writer::ResultWriters     m_writer;
};

} // namespace io::reader
} // namespace fem
