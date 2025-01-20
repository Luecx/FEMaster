#include "reader.h"
#include "../model/model.h"

namespace fem::reader {
void Reader::process_support() {
    // read SUPPORT, SUPPORT_COLLECTOR=xyz
    // NSET, lx, ly, lz
    // id, lx, ly, lz
    // ...
    // use empty field if no support
    auto supp_col    = m_current_line.require<std::string>("SUPPORT_COLLECTOR");
    auto orientation = m_current_line.parse<std::string>("ORIENTATION", "");

    m_model->_data->supp_cols.activate(supp_col);

    while (next_line().type() == DATA_LINE) {
        auto str = m_current_line.values()[0];
        StaticVector<6> constraint;

        for(Dim i = 0; i < 6; i++) {
            bool is_given = i < m_current_line.values().size() - 1;
            bool is_empty = !is_given || m_current_line.values()[i + 1].empty();

            if (is_given && !is_empty) {
                constraint(i) = std::stof(m_current_line.values()[i + 1]);
            } else {
                constraint(i) = NAN;
            }
        }

        if (m_model->_data->node_sets.has(str)) {
            m_model->add_support(str, constraint, orientation);
        } else {
            m_model->add_support(std::stoi(str), constraint, orientation);
        }
    }
}
} // namespace fem::reader
