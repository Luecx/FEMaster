//
// Created by f_eggers on 20.01.2025.
//
#include "reader.h"

#include "../model/model.h"

namespace fem::reader {
void Reader::process() {
    // read first line
    next_line();

    while (m_current_line.type() != END_OF_FILE) {
        while (m_current_line.type() != KEYWORD_LINE && m_current_line.type() != END_OF_FILE) {
            next_line();
        }

        if (m_current_line.type() == END_OF_FILE) {
            break;
        }

        logging::info(true, "Parsing: ", m_current_line.line());
        if (m_current_line.command() == "NODE") {
            process_nodes();
        } else if (m_current_line.command() == "ELEMENT") {
            process_elements();
        } else if (m_current_line.command() == "SURFACE") {
            process_surfaces();
        } else if (m_current_line.command() == "NSET") {
            process_nset();
        } else if (m_current_line.command() == "ELSET") {
            process_elset();
        } else if (m_current_line.command() == "SFSET") {
            process_sfset();
        } else if (m_current_line.command() == "MATERIAL") {
            process_material();
        } else if (m_current_line.command() == "ELASTIC") {
            process_elastic();
        } else if (m_current_line.command() == "DENSITY") {
            process_density();
        } else if (m_current_line.command() == "THERMALEXPANSION") {
            process_thermal_expansion();
        } else if (m_current_line.command() == "PROFILE") {
            process_profile();
        } else if (m_current_line.command() == "SOLIDSECTION") {
            process_solid_section();
        } else if (m_current_line.command() == "BEAMSECTION") {
            process_beam_section();
        } else if (m_current_line.command() == "SHELLSECTION") {
            process_shell_section();
        } else if (m_current_line.command() == "POINTMASSSECTION") {
            process_point_mass_section();
        } else if (m_current_line.command() == "SUPPORT") {
            process_support();
        } else if (m_current_line.command() == "TEMPERATURE") {
            process_temperature();
        } else if (m_current_line.command() == "ORIENTATION") {
            process_coordinate_system();
        } else if (m_current_line.command() == "CONNECTOR") {
            process_connector();
        } else if (m_current_line.command() == "COUPLING") {
            process_coupling();
        } else if (m_current_line.command() == "TIE") {
            process_tie();
        } else if (m_current_line.command() == "CLOAD") {
            process_cload();
        } else if (m_current_line.command() == "DLOAD") {
            process_dload();
        } else if (m_current_line.command() == "VLOAD") {
            process_vload();
        } else if (m_current_line.command() == "TLOAD") {
            process_tload();
        } else if (m_current_line.command() == "LOADCASE") {
            process_loadcase();
        } else if (m_current_line.command() == "EXIT") {
            std::exit(0); next_line();
        } else if (m_current_line.command() == "DEBUG") {
            std::cout << *m_model << std::endl; next_line();
        }

        else {
            logging::warning(false, "Unknown command: ", m_current_line.line());
            next_line();
        }
    }
}
} // namespace fem::reader
