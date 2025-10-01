#include "reader.h"

#include "../model/model.h"
#include "../constraints/connector.h"

namespace fem::reader {
void Reader::process_connector() {
    auto type = m_current_line.require<std::string>("TYPE");
    auto nset1 = m_current_line.require<std::string>("NSET1");
    auto nset2 = m_current_line.require<std::string>("NSET2");
    auto coord = m_current_line.require<std::string>("COORDINATESYSTEM");

    constraint::ConnectorType ctype = constraint::ConnectorType::None;
    if (type == "BEAM") {
        ctype = constraint::ConnectorType::Beam;
    } else if (type == "HINGE") {
        ctype = constraint::ConnectorType::Hinge;
    } else if (type == "CYLINDRICAL") {
        ctype = constraint::ConnectorType::Cylindrical;
    } else if (type == "TRANSLATOR") {
        ctype = constraint::ConnectorType::Translator;
    } else if (type == "JOIN") {
        ctype = constraint::ConnectorType::Join;
    } else {
        logging::error(false, "Unknown connector type: ", type);
    }

    m_model->add_connector(nset1, nset2, coord, ctype);

    next_line();
}
} // namespace fem::reader
