#include "reader.h"
#include "../model/model.h"

namespace fem::reader {
void Reader::process_solid_section() {
    auto mat = m_current_line.require<std::string>("MAT", "MATERIAL");
    auto els = m_current_line.require<std::string>("ELSET");
    m_model->solid_section(els, mat);
    next_line();
}
} // namespace fem::reader
