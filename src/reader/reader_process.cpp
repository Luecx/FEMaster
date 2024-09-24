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
        } else if (m_current_line.command() == "NSET") {
            process_nset();
        } else if (m_current_line.command() == "ELSET") {
            process_elset();
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
        } else if (m_current_line.command() == "COUPLING") {
            process_coupling();
        } else if (m_current_line.command() == "CLOAD") {
            process_cload();
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

void Reader::process_nset() {
    // read NSET, NAME=xyz
    auto setname = m_current_line.require<std::string>("NAME", "NSET");
    m_model->activate_node_set(setname);
    while (next_line().type() == DATA_LINE) {
        for (const auto& id : m_current_line.values()) {
            m_model->active_nodeset().push_back(std::stoi(id));
        }
    }
}

void Reader::process_elset() {
    // read ELSET, NAME=xyz
    m_model->activate_element_set(m_current_line.require<std::string>("NAME", "ELSET"));
    while (next_line().type() == DATA_LINE) {
        for (const auto& id : m_current_line.values()) {
            m_model->active_elemset().push_back(std::stoi(id));
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
            m_model->add_cload(str, StaticVector<3>(lx, ly, lz));
        } else {
            m_model->add_cload(std::stoi(m_current_line.values()[0]), StaticVector<3>(lx, ly, lz));
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
            m_model->add_vload(str, StaticVector<3>(lx, ly, lz));
        } else {
            m_model->add_vload(std::stoi(m_current_line.values()[0]), StaticVector<3>(lx, ly, lz));
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