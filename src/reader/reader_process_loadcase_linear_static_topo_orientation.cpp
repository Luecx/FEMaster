#include "reader.h"
#include "../loadcase/linear_static_topo.h"

namespace fem::reader {
void Reader::process_loadcase_linear_static_topo_orientation(fem::loadcase::LinearStaticTopo* lc) {
    while (next_line().type() == DATA_LINE) {
        auto id         = m_current_line.get_value(0, 0);
        auto angle_1    = m_current_line.get_value(1, 0.0f);
        auto angle_2    = m_current_line.get_value(2, 0.0f);
        auto angle_3    = m_current_line.get_value(3, 0.0f);
        lc->orientation(id, 0) = angle_1;
        lc->orientation(id, 1) = angle_2;
        lc->orientation(id, 2) = angle_3;
    }
}
} // namespace fem::reader
