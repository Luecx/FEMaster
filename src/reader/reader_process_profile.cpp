#include "reader.h"
#include "../model/model.h"

namespace fem::reader {
void Reader::process_profile() {
    auto mat = m_current_line.require<std::string>("NAME", "PROFILE");
    next_line();
    Precision A = m_current_line.get_value(0, 0.0f);
    Precision Iy = m_current_line.get_value(1, 0.0f);
    Precision Iz = m_current_line.get_value(2, 0.0f);
    Precision Jt = m_current_line.get_value(3, 0.0f);
    m_model->_data->profiles.activate(mat, A, Iy, Iz, Jt);
    next_line();
}
} // namespace fem::reader
