#include "reader.h"

#include "../model/model.h"
#include "../cos/rectangular_system.h"
#include "../cos/cylindrical_system.h"

namespace fem::reader {
void Reader::process_coordinate_system() {
    auto type       = m_current_line.require<std::string>("TYPE");
    auto definition = m_current_line.parse<std::string>("DEFINITION", "VECTOR");
    auto name       = m_current_line.require<std::string>("NAME");

    next_line();

    if (type == "RECTANGULAR") {
        if (m_current_line.count_values() == 3) {
            auto x = m_current_line.get_value(0, 0.0f);
            auto y = m_current_line.get_value(1, 0.0f);
            auto z = m_current_line.get_value(2, 0.0f);
            m_model->add_coordinate_system<cos::RectangularSystem>(name, Vec3{x,y,z});
        } else if (m_current_line.count_values() == 6) {
            auto x1 = m_current_line.get_value(0, 0.0f);
            auto y1 = m_current_line.get_value(1, 0.0f);
            auto z1 = m_current_line.get_value(2, 0.0f);
            auto x2 = m_current_line.get_value(3, 0.0f);
            auto y2 = m_current_line.get_value(4, 0.0f);
            auto z2 = m_current_line.get_value(5, 0.0f);
            m_model->add_coordinate_system<cos::RectangularSystem>(name, Vec3{x1,y1,z1}, Vec3{x2,y2,z2});
        } else if (m_current_line.count_values() == 9) {
            auto x1 = m_current_line.get_value(0, 0.0f);
            auto y1 = m_current_line.get_value(1, 0.0f);
            auto z1 = m_current_line.get_value(2, 0.0f);
            auto x2 = m_current_line.get_value(3, 0.0f);
            auto y2 = m_current_line.get_value(4, 0.0f);
            auto z2 = m_current_line.get_value(5, 0.0f);
            auto x3 = m_current_line.get_value(6, 0.0f);
            auto y3 = m_current_line.get_value(7, 0.0f);
            auto z3 = m_current_line.get_value(8, 0.0f);
            m_model->add_coordinate_system<cos::RectangularSystem>(name, Vec3{x1,y1,z1}, Vec3{x2,y2,z2}, Vec3{x3,y3,z3});
        } else {
            logging::error(false, "Cannot create coordinate system with", m_current_line.count_values(),
                " values. Needs to be 3, 6 or 9.");
        }
    } else if (type== "CYLINDRICAL") {
        if (m_current_line.count_values() == 9) {
            auto x1 = m_current_line.get_value(0, 0.0f);
            auto y1 = m_current_line.get_value(1, 0.0f);
            auto z1 = m_current_line.get_value(2, 0.0f);
            auto x2 = m_current_line.get_value(3, 0.0f);
            auto y2 = m_current_line.get_value(4, 0.0f);
            auto z2 = m_current_line.get_value(5, 0.0f);
            auto x3 = m_current_line.get_value(6, 0.0f);
            auto y3 = m_current_line.get_value(7, 0.0f);
            auto z3 = m_current_line.get_value(8, 0.0f);
            m_model->add_coordinate_system<cos::CylindricalSystem>(name, Vec3{x1,y1,z1}, Vec3{x2,y2,z2}, Vec3{x3,y3,z3});
        } else {
            logging::error(false, "Cannot create coordinate system with", m_current_line.count_values(),
                " values. Needs to be exactly 9.");
        }
    }

    else {
        logging::error(false, "Unknown coordinate system type: ", type);
    }
    next_line();
}
} // namespace fem::reader
