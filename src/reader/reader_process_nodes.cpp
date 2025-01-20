#include "reader.h"

#include "../model/model.h"

namespace fem::reader {
void Reader::process_nodes() {
    // read NODE
    // check if NSET is defined, if not, use NALL
    m_model->_data->node_sets.activate(m_current_line.parse<std::string>("NSET", "NALL"));

    while (next_line().type() == DATA_LINE) {
        int node_id = m_current_line.get_value(0, 0);
        Precision x = m_current_line.get_value(1, (Precision) 0);
        Precision y = m_current_line.get_value(2, (Precision) 0);
        Precision z = m_current_line.get_value(3, (Precision) 0);
        m_model->set_node(node_id, x, y, z);
    }
}
} // namespace fem::reader
