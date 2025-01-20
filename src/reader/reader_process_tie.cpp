#include "reader.h"
#include "../model/model.h"


namespace fem::reader {
void Reader::process_tie() {
    // read COUPLING, MASTER=xyz, SLAVE=xyz
    auto master_set = m_current_line.require<std::string>("MASTER");
    auto slave_set  = m_current_line.require<std::string>("SLAVE");
    auto adjust     = m_current_line.parse  <std::string>("ADJUST", "NO");
    auto distance   = m_current_line.require<Precision>("DISTANCE");

    m_model->add_tie(master_set, slave_set, distance, adjust == "YES");
    next_line();
}
} // namespace fem::reader
