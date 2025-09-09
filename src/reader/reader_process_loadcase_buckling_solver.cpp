#include "reader.h"
#include "../loadcase/linear_buckling.h"

namespace fem::reader {

void Reader::process_loadcase_linear_buckling_solver(fem::loadcase::LinearBuckling* lc) {
    // parse a single keyword/options line that may contain DEVICE and/or METHOD
    if (m_current_line.has_key("DEVICE")) {
        lc->device = m_current_line.parse<std::string>("DEVICE", "CPU") == "CPU"
                     ? solver::CPU
                     : solver::GPU;
    }
    if (m_current_line.has_key("METHOD")) {
        lc->method = m_current_line.parse<std::string>("METHOD", "DIRECT") == "DIRECT"
                     ? solver::DIRECT
                     : solver::INDIRECT;
    }
    next_line();
}

} // namespace fem::reader
