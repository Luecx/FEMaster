#include "reader.h"
#include "../loadcase/linear_static_topo.h"

namespace fem::reader {
void Reader::process_loadcase_linear_static_topo_density(fem::loadcase::LinearStaticTopo* lc) {
    while (next_line().type() == DATA_LINE) {
        auto id         = m_current_line.get_value(0, 0);
        auto ds         = m_current_line.get_value(1, 0.0f);

        lc->density(id) = ds;
    }
}
} // namespace fem::reader
