/**
* @file reader_loadcase_linear_buckling.cpp
 * @brief Reader extensions for LinearBuckling stiffness / geom file requests.
 *
 * @date    25.09.2025
 * @author  Finn
 */

#include "reader.h"
#include "../loadcase/linear_buckling.h"

namespace fem::reader {

void Reader::process_loadcase_linear_buckling_request_stiffness(fem::loadcase::LinearBuckling* lc) {
    if (m_current_line.has_key("FILE")) {
        lc->stiffness_file = m_current_line.parse<std::string>("FILE", "");
    } else {
        lc->stiffness_file = "buckling_stiffness_" + std::to_string(lc->get_id()) + ".mtx";
    }
    next_line();
}

void Reader::process_loadcase_linear_buckling_request_geom(fem::loadcase::LinearBuckling* lc) {
    if (m_current_line.has_key("FILE")) {
        lc->geom_file = m_current_line.parse<std::string>("FILE", "");
    } else {
        lc->geom_file = "buckling_geom_" + std::to_string(lc->get_id()) + ".mtx";
    }
    next_line();
}

} // namespace fem::reader
