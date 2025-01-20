#include "reader.h"
#include "../model/model.h"
#include "../material/material.h"

namespace fem::reader {
void Reader::process_density() {
    // read DENSITY
    next_line();
    auto v = m_current_line.get_value(0, 0.0f);
    m_model->_data->materials.get()->set_density(v);
    next_line();
}
} // namespace fem::reader
