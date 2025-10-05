#include "reader.h"

#include "../material/elasticity.h"
#include "../material/isotropic_elasticity.h"
#include "../material/orthotropic_elasticity.h"
#include "../model/model.h"

namespace fem::reader {
void Reader::process_elastic() {
    // read ELASTIC, TYPE=xyz
    auto type = m_current_line.require<std::string>("TYPE");
    next_line();
    if (type == "ISO" || type == "ISOTROPIC") {
        auto E = m_current_line.get_value(0, 0.0f);
        auto n = m_current_line.get_value(1, 0.0f);
        m_model->_data->materials.get()->set_elasticity<fem::material::IsotropicElasticity>(E, n);
    }
    if (type == "ORTHO" || type == "ORTHOTROPIC") {
        auto E1 = m_current_line.get_value(0, 0.0f);
        auto E2 = m_current_line.get_value(1, 0.0f);
        auto E3 = m_current_line.get_value(2, 0.0f);
        auto G23 = m_current_line.get_value(3, 0.0f);
        auto G13 = m_current_line.get_value(4, 0.0f);
        auto G12 = m_current_line.get_value(5, 0.0f);
        auto nu23 = m_current_line.get_value(6, 0.0f);
        auto nu13 = m_current_line.get_value(7, 0.0f);
        auto nu12 = m_current_line.get_value(8, 0.0f);
        m_model->_data->materials.get()->set_elasticity<fem::material::OrthotropicElasticity>(
            E1, E2, E3, G23, G13, G12, nu23, nu13, nu12);
    }
    next_line();
}
} // namespace fem::reader
