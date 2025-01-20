#include "reader.h"

#include "../model/model.h"
#include "../constraints/coupling.h"

namespace fem::reader {
void Reader::process_coupling() {
    // read COUPLING, MASTER=xyz, SLAVE=xyz, DOFS=xyz, TYPE=xyz
    auto master_set = m_current_line.require<std::string>("MASTER");
    auto type       = m_current_line.require<std::string>("TYPE");

    // either SLAVE or SURFACE can be defined. slave refers to a node set and surface to a surface set
    auto slave_set  = m_current_line.parse  <std::string>("SLAVE", "");
    auto surface    = m_current_line.parse  <std::string>("SFSET", "");

    next_line();

    Dofs dof_mask {false, false, false, false, false, false};
    for (Dim i = 0; i < 6; i++) {
        if (i < m_current_line.values().size()) {
            dof_mask(i) = std::stof(m_current_line.values()[i]) > 0;
        }
    }

    if (type == "KINEMATIC") {
        if (surface.empty()) {
            m_model->add_coupling(master_set, slave_set, dof_mask, constraint::CouplingType::KINEMATIC, false);
        } else {
            m_model->add_coupling(master_set, surface, dof_mask, constraint::CouplingType::KINEMATIC, true);
        }
    } else {
        logging::warning(false, "Unknown coupling type: ", type);
        logging::warning(false, "    Known coupling types: KINEMATIC");
    }
    next_line();
}
} // namespace fem::reader
