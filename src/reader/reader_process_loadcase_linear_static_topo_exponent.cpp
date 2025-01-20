#include "reader.h"
#include "../loadcase/linear_static_topo.h"

namespace fem::reader {
void Reader::process_loadcase_linear_static_topo_exponent(fem::loadcase::LinearStaticTopo* lc) {
    next_line();
    auto exp     = m_current_line.get_value(0, 0.0f);
    lc->exponent = exp;
    next_line();
}

}    // namespace fem::reader
