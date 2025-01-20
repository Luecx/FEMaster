#include "reader.h"

#include "../model/model.h"

namespace fem::reader {
void Reader::process_material() {
    // read MATERIAL, NAME=xyz
    m_model->_data->materials.activate(m_current_line.require<std::string>("NAME"));
    next_line();
}
} // namespace fem::reader
