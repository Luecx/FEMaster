#include "reader.h"
#include "../loadcase/linear_eigenfreq.h"

namespace fem::reader {
void Reader::process_loadcase_eigenfreq() {
    // read first line
    next_line();
    loadcase::LinearEigenfrequency lc {m_data.current_loadcase_num++, &m_writer, m_model.get(), 10};
    while (m_current_line.type() != END_OF_FILE) {
        while (m_current_line.type() != KEYWORD_LINE && m_current_line.type() != END_OF_FILE) {
            next_line();
        }

        logging::info(true, "   Parsing: ", m_current_line.line());
        if (m_current_line.command() == "NUMEIGENVALUES") {
            next_line();
            lc.num_eigenvalues = m_current_line.get_value(0, 10);
        } else if (m_current_line.command() == "SUPPORT") {
            process_loadcase_eigenfreq_support(&lc);
        } else if (m_current_line.command() == "END") {
            next_line();
            break;
        } else {
            logging::warning(false, "   Unknown command for eigenfrequency loadcase: ", m_current_line.line());
            next_line();
        }
    }

    lc.run();
}
} // namespace fem::reader
