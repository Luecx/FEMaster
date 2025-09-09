#include "reader.h"
#include "../loadcase/linear_buckling.h"

namespace fem::reader {

void Reader::process_loadcase_linear_buckling_support(fem::loadcase::LinearBuckling* lc) {
    // read subsequent DATA_LINEs as support set names
    while (next_line().type() == DATA_LINE) {
        for (const auto& str : m_current_line.values()) {
            lc->supps.push_back(str);
        }
    }
}

} // namespace fem::reader
