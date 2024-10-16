//
// Created by Finn Eggers on 05.09.23.
//
#include "../model/c3d10.h"
#include "../model/c3d15.h"
#include "../model/c3d20.h"
#include "../model/c3d20r.h"
#include "../model/c3d4.h"
#include "../model/c3d6.h"
#include "../model/c3d8.h"
#include "reader.h"
#include "../cos/rectangular_system.h"

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
        } else if (m_current_line.command() == "SOLIDSECTION") {
            process_solid_section();
        } else if (m_current_line.command() == "SUPPORT") {
            process_support();
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
        } else if (m_current_line.command() == "LOADCASE") {
            process_loadcase();
        } else {
            logging::warning(false, "Unknown command: ", m_current_line.line());
            next_line();
        }
    }
}

void Reader::process_nodes() {
    // read NODE
    // check if NSET is defined, if not, use NALL
    m_model->activate_node_set(m_current_line.parse<std::string>("NSET", ""));

    while (next_line().type() == DATA_LINE) {
        int node_id = std::stoi(m_current_line.values()[0]);
        Precision x = (Precision) std::stod(m_current_line.values()[1]);
        Precision y = (Precision) std::stod(m_current_line.values()[2]);
        Precision z = m_current_line.values().size() > 3 ? (Precision) std::stod(m_current_line.values()[3]) : 0;
        m_model->set_node(node_id, x, y, z);
    }
}

void Reader::process_elements() {
    // read ELEMENT
    // check if ELSET is defined, if not, use EALL
    m_model->activate_element_set(m_current_line.parse<std::string>("ELSET", ""));

    auto type = m_current_line.require<std::string>("TYPE");


    auto gather_values = [&](int req_values) -> std::vector<ID> {
        std::vector<ID> values;
        for(int i = 1; i < m_current_line.values().size(); i++){
            values.push_back(std::stoi(m_current_line.values()[i]));
        }
        while(values.size() < req_values){
            next_line();
            for(const auto& val : m_current_line.values()){
                values.push_back(std::stoi(val));
            }
        }
        return values;
    };

    while (next_line().type() == DATA_LINE) {
        int id = std::stoi(m_current_line.values()[0]);
        if (type == "C3D4") {
            auto values = gather_values(4);
            m_model->set_element<fem::model::C3D4>(id, values[0], values[1], values[2], values[3]);
        } else if (type == "C3D5") {
            auto values = gather_values(5);
            m_model->set_element<fem::model::C3D8>(id,
                                                   values[0], values[1], values[2], values[3],
                                                   values[4], values[4], values[4], values[4]);
        } else if (type == "C3D6") {
            auto values = gather_values(6);
            m_model->set_element<fem::model::C3D6>(id,
                                                   values[0], values[1], values[2], values[3],
                                                   values[4], values[5]);
        } else if (type == "C3D8") {
            auto values = gather_values(8);
            m_model->set_element<fem::model::C3D8>(id,
                                                   values[0], values[1], values[2], values[3],
                                                   values[4], values[5], values[6], values[7]);
        } else if (type == "C3D10") {
            auto values = gather_values(10);
            m_model->set_element<fem::model::C3D10>(id,
                                                    values[0], values[1], values[2], values[3],
                                                    values[4], values[5], values[6], values[7],
                                                    values[8], values[9]);
        } else if (type == "C3D15") {
            auto values = gather_values(15);
            m_model->set_element<fem::model::C3D15>(id,
                                                    values[0], values[1], values[2], values[3],
                                                    values[4], values[5], values[6], values[7],
                                                    values[8], values[9], values[10], values[11],
                                                    values[12], values[13], values[14]);
        } else if (type == "C3D20") {
            auto values = gather_values(20);
            m_model->set_element<fem::model::C3D20>(id,
                                                    values[0], values[1], values[2], values[3],
                                                    values[4], values[5], values[6], values[7],
                                                    values[8], values[9], values[10], values[11],
                                                    values[12], values[13], values[14], values[15],
                                                    values[16], values[17], values[18], values[19]);
        } else if (type == "C3D20R") {
            auto values = gather_values(20);
            m_model->set_element<fem::model::C3D20R>(id,
                                                    values[0], values[1], values[2], values[3],
                                                    values[4], values[5], values[6], values[7],
                                                    values[8], values[9], values[10], values[11],
                                                    values[12], values[13], values[14], values[15],
                                                    values[16], values[17], values[18], values[19]);
        } else {
            logging::warning(false, "Unknown element type ", type);
            return;
        }
    }
}

void Reader::process_surfaces() {
    // Activate the specified element set (if any) before processing surfaces.
    m_model->activate_surface_set(m_current_line.require<std::string>("SFSET", "NAME"));

    while(next_line().type() == DATA_LINE) {
        // Read the number of parts in the current line.
        int num_values = static_cast<int>(m_current_line.values().size());

        int id = -1; // Default value for ID if not specified.
        int elem_id = -1;
        int surf_side = -1;
        std::string set_name;

        if (num_values == 3) {
            // Three parts: ID, elementID, sideID
            id = std::stoi(m_current_line.values()[0]);
            elem_id = std::stoi(m_current_line.values()[1]);
            if (m_current_line.values()[2][0] == 'S') {
                surf_side = std::stoi(m_current_line.values()[2].substr(1));
            } else {
                surf_side = std::stoi(m_current_line.values()[2]);
            }

            m_model->set_surface(id, elem_id, surf_side);

        } else if (num_values == 2) {
            // Two parts: Either elementID/elementSet and sideID
            const std::string &first_value = m_current_line.values()[0];

            // Determine if the first value is a numeric element ID or a set name.
            bool is_numeric = std::all_of(first_value.begin(), first_value.end(), ::isdigit);

            if (is_numeric) {
                // First value is element ID
                elem_id = std::stoi(first_value);
                if (m_current_line.values()[1][0] == 'S') {
                    surf_side = std::stoi(m_current_line.values()[1].substr(1));
                } else {
                    surf_side = std::stoi(m_current_line.values()[1]);
                }

                m_model->set_surface(id, elem_id, surf_side);
            } else {
                // First value is element set name
                set_name = first_value;
                if (m_current_line.values()[1][0] == 'S') {
                    surf_side = std::stoi(m_current_line.values()[1].substr(1));
                } else {
                    surf_side = std::stoi(m_current_line.values()[1]);
                }

                // Use the set name to define the surface.
                m_model->set_surface(set_name, surf_side);
            }
        }
    }
}

void Reader::process_nset() {
    // read NSET, NAME=xyz
    bool generate = m_current_line.has_key("GENERATE");
    auto setname  = m_current_line.require<std::string>("NAME", "NSET");
    m_model->activate_node_set(setname);
    while (next_line().type() == DATA_LINE) {
        if (generate) {
            // require exactly 3 values, not more, not less
            if (m_current_line.count_values() != 3) {
                logging::error(false, "GENERATE requires exactly 3 values.");
            }
            ID id1 = std::stoi(m_current_line.values()[0]);
            ID id2 = std::stoi(m_current_line.values()[1]);
            ID inc = std::stoi(m_current_line.values()[2]);
            for (ID i = id1; i <= id2; i += inc) {
                m_model->active_nodeset().push_back(i);
            }
        } else {
            for (const auto& id : m_current_line.values()) {
                m_model->active_nodeset().push_back(std::stoi(id));
            }
        }
    }
}

void Reader::process_elset() {
    // read ELSET, NAME=xyz
    bool generate = m_current_line.has_key("GENERATE");

    m_model->activate_element_set(m_current_line.require<std::string>("NAME", "ELSET"));
    while (next_line().type() == DATA_LINE) {

        if (generate) {
            // require exactly 3 values, not more, not less
            if (m_current_line.count_values() != 3) {
                logging::error(false, "GENERATE requires exactly 3 values.");
            }
            ID id1 = std::stoi(m_current_line.values()[0]);
            ID id2 = std::stoi(m_current_line.values()[1]);
            ID inc = std::stoi(m_current_line.values()[2]);
            for (ID i = id1; i <= id2; i += inc) {
                m_model->active_elemset().push_back(i);
            }
        } else {
            for (const auto& id : m_current_line.values()) {
                m_model->active_elemset().push_back(std::stoi(id));
            }
        }
    }
}

void Reader::process_sfset() {
    // read ELSET, NAME=xyz
    bool generate = m_current_line.has_key("GENERATE");
    m_model->activate_surface_set(m_current_line.require<std::string>("NAME", "SFSET"));
    while (next_line().type() == DATA_LINE) {
        if (generate) {
            // require exactly 3 values, not more, not less
            if (m_current_line.count_values() != 3) {
                logging::error(false, "GENERATE requires exactly 3 values.");
            }
            ID id1 = std::stoi(m_current_line.values()[0]);
            ID id2 = std::stoi(m_current_line.values()[1]);
            ID inc = std::stoi(m_current_line.values()[2]);
            for (ID i = id1; i <= id2; i += inc) {
                m_model->active_surfset().push_back(i);
            }
        } else {
            for (const auto& id : m_current_line.values()) {
                m_model->active_surfset().push_back(std::stoi(id));
            }
        }
    }
}

void Reader::process_material() {
    // read MATERIAL, NAME=xyz
    m_model->activate_material(m_current_line.require<std::string>("NAME"));
    next_line();
}

void Reader::process_elastic() {
    // read ELASTIC, TYPE=xyz
    auto type = m_current_line.require<std::string>("TYPE");
    next_line();
    if (type == "ISO" || type == "ISOTROPIC") {
        auto E = std::stof(m_current_line.values()[0]);
        auto n = std::stof(m_current_line.values()[1]);
        m_model->active_material().set_elasticity<fem::material::IsotropicElasticity>(E, n);
    }
    if (type == "ORTHO" || type == "ORTHOTROPIC") {
        auto E1 = std::stof(m_current_line.values()[0]);
        auto E2 = std::stof(m_current_line.values()[1]);
        auto E3 = std::stof(m_current_line.values()[2]);
        auto G23 = std::stof(m_current_line.values()[3]);
        auto G13 = std::stof(m_current_line.values()[4]);
        auto G12 = std::stof(m_current_line.values()[5]);
        auto nu23 = std::stof(m_current_line.values()[6]);
        auto nu13 = std::stof(m_current_line.values()[7]);
        auto nu12 = std::stof(m_current_line.values()[8]);
        m_model->active_material().set_elasticity<fem::material::OrthotropicElasticity>(E1, E2, E3, G23, G13, G12, nu23, nu13, nu12);
    }
    next_line();
}

void Reader::process_density() {
    // read DENSITY
    next_line();
    auto v = std::stof(m_current_line.values()[0]);
    m_model->active_material().set_density(v);
    next_line();
}

void Reader::process_solid_section() {
    auto mat = m_current_line.require<std::string>("MAT", "MATERIAL");
    auto els = m_current_line.require<std::string>("ELSET");
    m_model->solid_section(els, mat);
    next_line();
}

void Reader::process_cload() {
    // read CLOAD, LOAD_COLLECTOR=xyz
    // NSET, lx, ly, lz
    // id, lx, ly, lz
    // ...
    m_model->activate_load_set(m_current_line.require<std::string>("LOAD_COLLECTOR"));
    while (next_line().type() == DATA_LINE) {
        auto str = m_current_line.values()[0];
        auto lx  = std::stof(m_current_line.values()[1]);
        auto ly  = std::stof(m_current_line.values()[2]);
        auto lz  = m_current_line.values().size() > 3 ? std::stof(m_current_line.values()[3]) : 0;

        if (m_model->nodesets().has(str)) {
            m_model->add_cload(str, Vec3(lx, ly, lz));
        } else {
            m_model->add_cload(std::stoi(m_current_line.values()[0]), Vec3(lx, ly, lz));
        }
    }
}

void Reader::process_dload() {
    // read DLOAD, LOAD_COLLECTOR=xyz
    // SFSET, lx, ly, lz
    // id, lx, ly, lz
    // ...
    m_model->activate_load_set(m_current_line.require<std::string>("LOAD_COLLECTOR"));
    while (next_line().type() == DATA_LINE) {
        auto str = m_current_line.values()[0];
        auto lx  = std::stof(m_current_line.values()[1]);
        auto ly  = std::stof(m_current_line.values()[2]);
        auto lz  = m_current_line.values().size() > 3 ? std::stof(m_current_line.values()[3]) : 0;

        if (m_model->surfsets().has(str)) {
            m_model->add_dload(str, Vec3(lx, ly, lz));
        } else {
            m_model->add_dload(std::stoi(m_current_line.values()[0]), Vec3(lx, ly, lz));
        }
    }
}

void Reader::process_vload() {
    // read VLOAD, LOAD_COLLECTOR=xyz
    // NSET, lx, ly, lz
    // id, lx, ly, lz
    // ...
    m_model->activate_load_set(m_current_line.require<std::string>("LOAD_COLLECTOR"));
    while (next_line().type() == DATA_LINE) {
        auto str = m_current_line.values()[0];
        auto lx  = std::stof(m_current_line.values()[1]);
        auto ly  = std::stof(m_current_line.values()[2]);
        auto lz  = m_current_line.values().size() > 3 ? std::stof(m_current_line.values()[3]) : 0;

        if (m_model->elemsets().has(str)) {
            m_model->add_vload(str, Vec3(lx, ly, lz));
        } else {
            m_model->add_vload(std::stoi(m_current_line.values()[0]), Vec3(lx, ly, lz));
        }
    }
}

void Reader::process_support() {
    // read SUPPORT, SUPPORT_COLLECTOR=xyz
    // NSET, lx, ly, lz
    // id, lx, ly, lz
    // ...
    // use empty field if no support
    m_model->activate_support_set(m_current_line.require<std::string>("SUPPORT_COLLECTOR"));

    while (next_line().type() == DATA_LINE) {
        auto str = m_current_line.values()[0];
        StaticVector<6> constraint;

        for(int i = 0; i < 6; i++) {
            bool is_given = i < m_current_line.values().size() - 1;
            bool is_empty = !is_given || m_current_line.values()[i + 1].empty();

            if (is_given && !is_empty) {
                constraint(i) = std::stof(m_current_line.values()[i + 1]);
            } else {
                constraint(i) = NAN;
            }
        }

        if (m_model->nodesets().has(str)) {
            m_model->add_support(str, constraint);
        } else {
            m_model->add_support(std::stoi(str), constraint);
        }
    }
}

void Reader::process_coordinate_system() {
    auto type       = m_current_line.require<std::string>("TYPE");
    auto definition = m_current_line.parse<std::string>("DEFINITION", "VECTOR");
    auto name       = m_current_line.require<std::string>("NAME");

    next_line();

    if (type == "RECTANGULAR") {
        if (m_current_line.count_values() == 3) {
            auto x = std::stof(m_current_line.values()[0]);
            auto y = std::stof(m_current_line.values()[1]);
            auto z = std::stof(m_current_line.values()[2]);
            m_model->add_coordinate_system<cos::RectangularSystem>(name, Vec3{x,y,z});
        } else if (m_current_line.count_values() == 6) {
            auto x1 = std::stof(m_current_line.values()[0]);
            auto y1 = std::stof(m_current_line.values()[1]);
            auto z1 = std::stof(m_current_line.values()[2]);
            auto x2 = std::stof(m_current_line.values()[3]);
            auto y2 = std::stof(m_current_line.values()[4]);
            auto z2 = std::stof(m_current_line.values()[5]);
            m_model->add_coordinate_system<cos::RectangularSystem>(name, Vec3{x1,y1,z1}, Vec3{x2,y2,z2});
        } else if (m_current_line.count_values() == 6) {
            auto x1 = std::stof(m_current_line.values()[0]);
            auto y1 = std::stof(m_current_line.values()[1]);
            auto z1 = std::stof(m_current_line.values()[2]);
            auto x2 = std::stof(m_current_line.values()[3]);
            auto y2 = std::stof(m_current_line.values()[4]);
            auto z2 = std::stof(m_current_line.values()[5]);
            auto x3 = std::stof(m_current_line.values()[6]);
            auto y3 = std::stof(m_current_line.values()[7]);
            auto z3 = std::stof(m_current_line.values()[8]);
            m_model->add_coordinate_system<cos::RectangularSystem>(name, Vec3{x1,y1,z1}, Vec3{x2,y2,z2}, Vec3{x3,y3,z3});
        } else {
            logging::error(false, "Cannot create coordinate system with", m_current_line.count_values(),
                " values. Needs to be 3, 6 or 9.");
        }
    } else {
        logging::error(false, "Unknown coordinate system type: ", type);
    }
    next_line();
}

void Reader::process_connector() {
    auto type = m_current_line.require<std::string>("TYPE");
    auto nset1 = m_current_line.require<std::string>("NSET1");
    auto nset2 = m_current_line.require<std::string>("NSET2");
    auto coord = m_current_line.require<std::string>("COORDINATESYSTEM");

    ConnectorType ctype = ConnectorType::None;
    if (type == "BEAM") {
        ctype = ConnectorType::Beam;
    } else if (type == "HINGE") {
        ctype = ConnectorType::Hinge;
    } else if (type == "CYLINDRICAL") {
        ctype = ConnectorType::Cylindrical;
    } else if (type == "TRANSLATOR") {
        ctype = ConnectorType::Translator;
    } else {
        logging::error(false, "Unknown connector type: ", type);
    }

    m_model->add_connector(nset1, nset2, coord, ctype);

    next_line();
}


void Reader::process_coupling() {
    // read COUPLING, MASTER=xyz, SLAVE=xyz, DOFS=xyz, TYPE=xyz
    auto master_set = m_current_line.require<std::string>("MASTER");
    auto slave_set  = m_current_line.require<std::string>("SLAVE");
    auto type       = m_current_line.require<std::string>("TYPE");

    next_line();

    Dofs dof_mask {false, false, false, false, false, false};
    for (int i = 0; i < 6; i++) {
        if (i < m_current_line.values().size()) {
            dof_mask(i) = std::stof(m_current_line.values()[i]) > 0;
        }
    }

    if (type == "KINEMATIC") {
        m_model->add_coupling(master_set, slave_set, dof_mask, CouplingType::KINEMATIC);
    } else {
        logging::warning(false, "Unknown coupling type: ", type);
        logging::warning(false, "    Known coupling types: KINEMATIC");
    }
    next_line();
}

void Reader::process_tie() {
    // read COUPLING, MASTER=xyz, SLAVE=xyz
    auto master_set = m_current_line.require<std::string>("MASTER");
    auto slave_set  = m_current_line.require<std::string>("SLAVE");
    auto adjust     = m_current_line.parse  <std::string>("ADJUST", "NO");
    auto distance   = m_current_line.require<Precision>("DISTANCE");

    m_model->add_tie(master_set, slave_set, distance, adjust == "YES");
    next_line();
}

void Reader::process_loadcase() {
    if (m_current_line.require<std::string>("TYPE") == "LINEARSTATIC") {
        process_loadcase_linear_static();
    } else if (m_current_line.require<std::string>("TYPE") == "LINEARSTATICTOPO") {
        process_loadcase_linear_static_topo();
    } else if (m_current_line.require<std::string>("TYPE") == "EIGENFREQ") {
        process_loadcase_eigenfreq();
    } else {
        logging::warning(false, "unknown loadcase type: ", m_current_line.require<std::string>("TYPE"));
    }
}

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

void Reader::process_loadcase_linear_static_topo() {
    // read first line
    next_line();
    loadcase::LinearStaticTopo lc {m_data.current_loadcase_num++, &m_writer, m_model.get()};
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
        } else if (m_current_line.command() == "EXPONENT") {
            process_loadcase_linear_static_topo_exponent(&lc);
        } else if (m_current_line.command() == "DENSITY") {
            process_loadcase_linear_static_topo_density(&lc);
        } else if (m_current_line.command() == "END") {
            next_line();
            break;
        } else {
            logging::warning(false, "   Unknown command for linear static topology loadcase: ", m_current_line.line());
            next_line();
        }
    }

    lc.run();
}

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

void Reader::process_loadcase_eigenfreq_support(fem::loadcase::LinearEigenfrequency* lc) {
    while (next_line().type() == DATA_LINE) {
        for (auto str : m_current_line.values()) {
            lc->supps.push_back(str);
        }
    }
}

void Reader::process_loadcase_linear_static_support(fem::loadcase::LinearStatic* lc) {
    while (next_line().type() == DATA_LINE) {
        for (auto str : m_current_line.values()) {
            lc->supps.push_back(str);
        }
    }
}
void Reader::process_loadcase_linear_static_load(fem::loadcase::LinearStatic* lc) {
    while (next_line().type() == DATA_LINE) {
        for (auto str : m_current_line.values()) {
            lc->loads.push_back(str);
        }
    }
}
void Reader::process_loadcase_linear_static_solver(fem::loadcase::LinearStatic* lc) {
    if (m_current_line.has_key("DEVICE")) {
        lc->device = m_current_line.parse<std::string>("DEVICE", "CPU") == "CPU" ? solver::CPU : solver::GPU;
    }
    if (m_current_line.has_key("METHOD")) {
        lc->method =
            m_current_line.parse<std::string>("METHOD", "DIRECT") == "DIRECT" ? solver::DIRECT : solver::INDIRECT;
    }
    next_line();
}

void Reader::process_loadcase_linear_static_request_stiffness(fem::loadcase::LinearStatic* lc) {
    if (m_current_line.has_key("FILE")) {
        lc->stiffness_file = m_current_line.parse<std::string>("FILE", "");
    } else {
        lc->stiffness_file = "stiffness_" + std::to_string(lc->id()) + ".txt";
    }
    next_line();
}

void Reader::process_loadcase_linear_static_topo_density(fem::loadcase::LinearStaticTopo* lc) {
    while (next_line().type() == DATA_LINE) {
        auto id         = std::stoi(m_current_line.values()[0]);
        auto ds         = std::stof(m_current_line.values()[1]);

        lc->density(id) = ds;
    }
}
void Reader::process_loadcase_linear_static_topo_exponent(fem::loadcase::LinearStaticTopo* lc) {
    next_line();
    auto exp     = std::stof(m_current_line.values()[0]);
    lc->exponent = exp;
    next_line();
}

}    // namespace fem::reader