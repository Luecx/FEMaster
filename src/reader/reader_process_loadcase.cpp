#include "reader.h"

namespace fem::reader {
void Reader::process_loadcase() {
    if (m_current_line.require<std::string>("TYPE") == "LINEARSTATIC") {
        process_loadcase_linear_static();
    } else if (m_current_line.require<std::string>("TYPE") == "LINEARSTATICTOPO") {
        process_loadcase_linear_static_topo();
    } else if (m_current_line.require<std::string>("TYPE") == "EIGENFREQ") {
        process_loadcase_eigenfreq();
    } else {
        logging::warning(false, "unknown loadcase type: ", m_current_line.require<std::string>("TYPE"));
    }
}
} // namespace fem::reader
