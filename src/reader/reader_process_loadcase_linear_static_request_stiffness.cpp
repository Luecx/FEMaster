#include "reader.h"
#include "../loadcase/linear_static.h"

namespace fem::reader {
void Reader::process_loadcase_linear_static_request_stiffness(fem::loadcase::LinearStatic* lc) {
    if (m_current_line.has_key("FILE")) {
        lc->stiffness_file = m_current_line.parse<std::string>("FILE", "");
    } else {
        lc->stiffness_file = "stiffness_" + std::to_string(lc->get_id()) + ".txt";
    }
    next_line();
}
} // namespace fem::reader
