#include "reader.h"
#include "../model/model.h"

namespace fem::reader {
void Reader::process_shell_section() {
    auto mat = m_current_line.require<std::string>("MAT", "MATERIAL");
    auto els = m_current_line.require<std::string>("ELSET");
    next_line();
    auto thi = m_current_line.get_value(0, 1.0f);
    m_model->shell_section(els, mat, thi);
    next_line();
}
} // namespace fem::reader
