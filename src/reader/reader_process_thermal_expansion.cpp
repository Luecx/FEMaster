#include "reader.h"
#include "../model/model.h"

namespace fem::reader {
void Reader::process_thermal_expansion() {
    // read THERMAL EXPANSION
    next_line();
    auto v = m_current_line.get_value(0, 0.0f);
    m_model->_data->materials.get()->set_thermal_expansion(v);
    next_line();
}
} // namespace fem::reader
