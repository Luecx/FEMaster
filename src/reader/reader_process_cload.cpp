#include "reader.h"

#include "../model/model.h"

namespace fem::reader {
void Reader::process_cload() {
    // read CLOAD, LOAD_COLLECTOR=xyz
    // NSET, lx, ly, lz
    // id, lx, ly, lz
    // ...
    m_model->_data->load_cols.activate(m_current_line.require<std::string>("LOAD_COLLECTOR"));
    while (next_line().type() == DATA_LINE) {
        auto str = m_current_line.values()[0];
        auto lx  = m_current_line.get_value(1, 0.0f);
        auto ly  = m_current_line.get_value(2, 0.0f);
        auto lz  = m_current_line.get_value(3, 0.0f);

        auto mx  = m_current_line.get_value(4, 0.0f);
        auto my  = m_current_line.get_value(5, 0.0f);
        auto mz  = m_current_line.get_value(6, 0.0f);

        if (m_model->_data->node_sets.has(str)) {
            m_model->add_cload(str, Vec6(lx, ly, lz, mx, my, mz));
        } else {
            m_model->add_cload(std::stoi(m_current_line.values()[0]), Vec6(lx, ly, lz, mx, my, mz));
        }
    }
}
} // namespace fem::reader
