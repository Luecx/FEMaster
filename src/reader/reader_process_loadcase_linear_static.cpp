#include "reader.h"
#include "../loadcase/linear_static.h"

namespace fem::reader {
void Reader::process_loadcase_linear_static() {
    // read first line
    next_line();
    loadcase::LinearStatic lc {m_data.current_loadcase_num++, &m_writer, m_model.get()};
    while (m_current_line.type() != END_OF_FILE) {
        while (m_current_line.type() != KEYWORD_LINE && m_current_line.type() != END_OF_FILE) {
            next_line();
        }

        logging::info(true, "   Parsing: ", m_current_line.line());
        if (m_current_line.command() == "SUPPORT") {
            process_loadcase_linear_static_support(&lc);
        } else if (m_current_line.command() == "LOAD") {
            process_loadcase_linear_static_load(&lc);
        } else if (m_current_line.command() == "SOLVER") {
            process_loadcase_linear_static_solver(&lc);
        } else if (m_current_line.command() == "REQUESTSTIFFNESS") {
            process_loadcase_linear_static_request_stiffness(&lc);
        }else if (m_current_line.command() == "END") {
            next_line();
            break;
        } else {
            logging::warning(false, "   Unknown command for linear static loadcase: ", m_current_line.line());
            next_line();
        }
    }

    lc.run();
}
} // namespace fem::reader
