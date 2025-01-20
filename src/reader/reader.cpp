//
// Created by Luecx on 06.09.2023.
//
#include "reader.h"
#include "../model/model.h"

namespace fem::reader{

Line& Reader::next_line() {
    m_current_line = m_file.next_line();
    return m_current_line;
}

Reader::Reader(const std::string& file_path, const std::string& output_path)
    :
    m_file(file_path),
    m_file_path(file_path),
    m_writer(output_path),
    m_output_path(output_path){}

void Reader::read() {
    m_file   = File(m_file_path);
    m_writer = Writer(m_output_path);

    // First stage
    analyse();

    // create the model
    m_model = std::make_shared<fem::model::Model>(m_data.highest_node_id + 1,
                                                  m_data.highest_element_id + 1,
                                                  m_data.highest_surface_id + 1);

    // Reset File for the next stage
    m_file = File(m_file_path);

    // Second stage
    process();

    // close the writer
    m_writer.close();
}


}
