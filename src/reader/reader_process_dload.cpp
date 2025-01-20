#include "reader.h"
#include "../model/model.h"

namespace fem::reader {
void Reader::process_dload() {
    // read DLOAD, LOAD_COLLECTOR=xyz
    // SFSET, lx, ly, lz
    // id, lx, ly, lz
    // ...
    m_model->_data->load_cols.activate(m_current_line.require<std::string>("LOAD_COLLECTOR"));

    while (next_line().type() == DATA_LINE) {
        auto str = m_current_line.values()[0];
        auto lx  = m_current_line.get_value(1, 0.0f);
        auto ly  = m_current_line.get_value(2, 0.0f);
        auto lz  = m_current_line.values().size() > 3 ? std::stof(m_current_line.values()[3]) : 0;

        if (m_model->_data->surface_sets.has(str)) {
            m_model->add_dload(str, Vec3(lx, ly, lz));
        } else {
            m_model->add_dload(std::stoi(m_current_line.values()[0]), Vec3(lx, ly, lz));
        }
    }
}
} // namespace fem::reader
