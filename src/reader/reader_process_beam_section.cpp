#include "reader.h"

#include "../model/model.h"

namespace fem::reader {
void Reader::process_beam_section() {
    auto mat = m_current_line.require<std::string>("MAT", "MATERIAL");
    auto els = m_current_line.require<std::string>("ELSET");
    auto profile = m_current_line.require<std::string>("PROFILE");
    next_line();
    Precision n1x = m_current_line.get_value(0, 0.0f);
    Precision n1y = m_current_line.get_value(1, 0.0f);
    Precision n1z = m_current_line.get_value(2, 0.0f);
    Vec3 n1 = Vec3(n1x, n1y, n1z);
    m_model->beam_section(els, mat, profile, n1);
    next_line();
}
} // namespace fem::reader
