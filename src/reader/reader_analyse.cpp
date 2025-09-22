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
    // Determine the required number of nodes for this element type
    auto type = m_current_line.require<std::string>("TYPE");

    while (next_line().type() == DATA_LINE) {
        // Get the element ID from the current line
        int element_id = std::stoi(m_current_line.values()[0]);
        m_data.highest_element_id = std::max(m_data.highest_element_id, element_id);

        Index req_values;
        if (type == "C3D4") {
            req_values = 4;
        } else if (type == "C3D5") {
            req_values = 5;
        } else if (type == "C3D6") {
            req_values = 6;
        } else if (type == "C3D8") {
            req_values = 8;
        } else if (type == "C3D10") {
            req_values = 10;
        } else if (type == "C3D15") {
            req_values = 15;
        } else if (type == "C3D20" || type == "C3D20R") {
            req_values = 20;
        } else if (type == "B33") {
            req_values = 2;
        } else if (type == "T3") {
            req_values = 2;
        } else if (type == "S3") {
            req_values = 3;
        } else if (type == "S4") {
            req_values = 4;
        } else if (type == "MITC4") {
            req_values = 4;
        } else if (type == "S6") {
            req_values = 6;
        } else if (type == "S8") {
            req_values = 8;
        } else if (type == "P") {
            req_values = 1;
        }else {
            logging::warning(false, "Unknown element type ", type);
            return;
        }

        // Gather values, including across multiple lines if necessary
        std::vector<ID> values;
        for (Index i = 1; i < m_current_line.values().size(); i++) {
            values.push_back(std::stoi(m_current_line.values()[i]));
        }
        while (values.size() < req_values) {
            next_line();
            for (const auto& val : m_current_line.values()) {
                values.push_back(std::stoi(val));
            }
        }
    }
}

void Reader::analyse_surfaces() {
    while (next_line().type() == DATA_LINE) {
        try {
            int surface_id            = std::stoi(m_current_line.values()[0]);
            m_data.highest_surface_id = std::max(m_data.highest_surface_id, surface_id);
        } catch (const std::invalid_argument& e) {
        }
    }
}


}    // namespace fem::reader