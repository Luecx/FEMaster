#include "reader.h"
#include "../model/model.h"

namespace fem::reader {
void Reader::process_tload() {
    // read TLOAD, LOAD_COLLECTOR=xyz, temperature field=xyz, reference temperature=xyz
    auto lod_col  = m_current_line.require<std::string>("LOAD_COLLECTOR");
    auto temp_fi  = m_current_line.require<std::string>("TEMPERATUREFIELD");
    auto ref_temp = m_current_line.require<Precision  >("REFERENCETEMPERATURE");

    m_model->_data->load_cols.activate(lod_col);
    m_model->add_tload(temp_fi, ref_temp);
    next_line();
}
} // namespace fem::reader
