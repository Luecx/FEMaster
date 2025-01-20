#include "reader.h"
#include "../model/model.h"

namespace fem::reader {
void Reader::process_sfset() {
    // read ELSET, NAME=xyz
    bool generate = m_current_line.has_key("GENERATE");
    m_model->_data->surface_sets.activate(m_current_line.require<std::string>("NAME", "SFSET"));

    while (next_line().type() == DATA_LINE) {
        if (generate) {
            // require exactly 2 or 3 values, not more, not less
            logging::warning(m_current_line.count_values() == 2
                            || m_current_line.count_values() == 3, "GENERATE requires exactly 2 or 3 values.");

            ID id1 = m_current_line.get_value(0, 0);
            ID id2 = m_current_line.get_value(1, 0);
            ID inc = m_current_line.get_value(2, 1);
            for (ID i = id1; i <= id2; i += inc) {
                m_model->_data->surface_sets.get()->add(i);
            }
        } else {
            for (const auto& id : m_current_line.values()) {
                m_model->_data->surface_sets.get()->add(std::stoi(id));
            }
        }
    }
}
} // namespace fem::reader
