//
// Created by Finn Eggers on 05.09.23.
//
#include "reader.h"

namespace fem::reader {

void Reader::analyse() {
    // read first line
    next_line();

    while (m_current_line.type() != END_OF_FILE) {
        while (m_current_line.type() != KEYWORD_LINE && m_current_line.type() != END_OF_FILE) {
            next_line();
        }

        if (m_current_line.command() == "NODE") {
            analyse_nodes();
        } else if (m_current_line.command() == "ELEMENT") {
            analyse_elements();
        } else if (m_current_line.command() == "SURFACE") {
            analyse_surfaces();
        } else {
            next_line();
        }
    }
}
void Reader::analyse_nodes() {
    while (next_line().type() == DATA_LINE) {
        int node_id            = std::stoi(m_current_line.values()[0]);
        m_data.highest_node_id = std::max(m_data.highest_node_id, node_id);
    }
}
void Reader::analyse_elements() {
    while (next_line().type() == DATA_LINE) {
        int element_id            = std::stoi(m_current_line.values()[0]);
        m_data.highest_element_id = std::max(m_data.highest_element_id, element_id);
    }
}
void Reader::analyse_surfaces() {
    while (next_line().type() == DATA_LINE) {
        int surface_id            = std::stoi(m_current_line.values()[0]);
        m_data.highest_surface_id = std::max(m_data.highest_surface_id, surface_id);
    }
}


}    // namespace fem::reader