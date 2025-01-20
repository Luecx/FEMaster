#include "reader.h"

#include "../section/section_point_mass.h"
#include "../model/model.h"


namespace fem::reader {
void Reader::process_point_mass_section() {
    auto els = m_current_line.require<std::string>("ELSET");
    next_line();
    Precision mass = 0;
    Vec3 inertia = Vec3::Zero();
    Vec3 spring  = Vec3::Zero();
    Vec3 rotary_spring = Vec3::Zero();

    if (m_current_line.type() == DATA_LINE) {
        mass = m_current_line.get_value(0, 0.0f);
        next_line();
    }

    if (m_current_line.type() == DATA_LINE) {
        inertia(0) = m_current_line.get_value(0, 0.0f);
        inertia(1) = m_current_line.get_value(1, 0.0f);
        inertia(2) = m_current_line.get_value(2, 0.0f);
        next_line();
    }

    if (m_current_line.type() == DATA_LINE) {
        spring(0) = m_current_line.get_value(0, 0.0f);
        spring(1) = m_current_line.get_value(1, 0.0f);
        spring(2) = m_current_line.get_value(2, 0.0f);
        next_line();
    }

    if (m_current_line.type() == DATA_LINE) {
        rotary_spring(0) = m_current_line.get_value(0, 0.0f);
        rotary_spring(1) = m_current_line.get_value(1, 0.0f);
        rotary_spring(2) = m_current_line.get_value(2, 0.0f);
        next_line();
    }

    m_model->point_mass_section(els, mass, inertia, spring, rotary_spring);
}
} // namespace fem::reader
