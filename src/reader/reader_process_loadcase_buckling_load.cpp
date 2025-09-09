#include "reader.h"
#include "../loadcase/linear_buckling.h"

namespace fem::reader {

void Reader::process_loadcase_linear_buckling_load(fem::loadcase::LinearBuckling* lc) {
    // read subsequent DATA_LINEs as load set names
    while (next_line().type() == DATA_LINE) {
        for (const auto& str : m_current_line.values()) {
            lc->loads.push_back(str);
        }
    }
}

} // namespace fem::reader
