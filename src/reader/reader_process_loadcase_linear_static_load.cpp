#include "reader.h"
#include "../loadcase/linear_static.h"

namespace fem::reader {
void Reader::process_loadcase_linear_static_load(fem::loadcase::LinearStatic* lc) {
    while (next_line().type() == DATA_LINE) {
        for (auto str : m_current_line.values()) {
            lc->loads.push_back(str);
        }
    }
}
} // namespace fem::reader
