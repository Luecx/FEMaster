#include "reader.h"
#include "../model/model.h"

namespace fem::reader {
void Reader::process_surfaces() {
    // Activate the specified element set (if any) before processing surfaces.
    m_model->_data->surface_sets.activate(m_current_line.require<std::string>("SFSET", "NAME"));

    while(next_line().type() == DATA_LINE) {
        // Read the number of parts in the current line.
        int num_values = static_cast<int>(m_current_line.values().size());

        int id = -1; // Default value for ID if not specified.
        int elem_id = -1;
        int surf_side = -1;
        std::string set_name;

        if (num_values == 3) {
            // Three parts: ID, elementID, sideID
            id      = m_current_line.get_value(0, 0);
            elem_id = m_current_line.get_value(1, 0);
            if (m_current_line.values()[2][0] == 'S') {
                surf_side = std::stoi(m_current_line.values()[2].substr(1));
            } else {
                surf_side = m_current_line.get_value(2, 0);
            }

            m_model->set_surface(id, elem_id, surf_side);

        } else if (num_values == 2) {
            // Two parts: Either elementID/elementSet and sideID
            const std::string &first_value = m_current_line.values()[0];

            // Determine if the first value is a numeric element ID or a set name.
            bool is_numeric = std::all_of(first_value.begin(), first_value.end(), ::isdigit);

            if (is_numeric) {
                // First value is element ID
                elem_id = std::stoi(first_value);
                if (m_current_line.values()[1][0] == 'S') {
                    surf_side = std::stoi(m_current_line.values()[1].substr(1));
                } else {
                    surf_side = m_current_line.get_value(1, 0);
                }

                m_model->set_surface(id, elem_id, surf_side);
            } else {
                // First value is element set name
                set_name = first_value;
                if (m_current_line.values()[1][0] == 'S') {
                    surf_side = std::stoi(m_current_line.values()[1].substr(1));
                } else {
                    surf_side = m_current_line.get_value(1, 0);
                }

                // Use the set name to define the surface.
                m_model->set_surface(set_name, surf_side);
            }
        }
    }
}
} // namespace fem::reader
