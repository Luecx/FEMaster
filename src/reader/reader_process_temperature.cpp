#include "reader.h"
#include "../model/model.h"

namespace fem::reader {
void Reader::process_temperature() {
    // read TEMPERATURE, NAME=xyz
    // NSET, value
    // id, value
    // ...
    auto name = m_current_line.require<std::string>("NAME");

    while (next_line().type() == DATA_LINE) {
        auto str   = m_current_line.values()[0];
        auto value = m_current_line.get_value(1, 0.0f);

        if (m_model->_data->node_sets.has(str)) {
            for (auto id : *m_model->_data->node_sets.get(str)) {
                m_model->set_field_temperature(name, id, value);
            }
        } else {
            m_model->set_field_temperature(name, std::stoi(str), value);
        }
    }
}
} // namespace fem::reader
