#include "reader.h"
#include "../model/model.h"

namespace fem::reader {
void Reader::process_pload() {
    // read DLOAD, LOAD_COLLECTOR=xyz
    // SFSET, lx, ly, lz
    // id, lx, ly, lz
    // ...
    m_model->_data->load_cols.activate(m_current_line.require<std::string>("LOAD_COLLECTOR"));

    while (next_line().type() == DATA_LINE) {
        auto str = m_current_line.values()[0];
        auto lv  = m_current_line.get_value(1, 0.0f);

        if (m_model->_data->surface_sets.has(str)) {
            m_model->add_pload(str, lv);
        } else {
            m_model->add_pload(std::stoi(m_current_line.values()[0]), lv);
        }
    }
}
} // namespace fem::reader
