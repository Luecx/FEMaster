#include "reader.h"
#include "../loadcase/linear_eigenfreq.h"

namespace fem::reader {
void Reader::process_loadcase_eigenfreq_support(fem::loadcase::LinearEigenfrequency* lc) {
    while (next_line().type() == DATA_LINE) {
        for (auto str : m_current_line.values()) {
            lc->supps.push_back(str);
        }
    }
}
} // namespace fem::reader
