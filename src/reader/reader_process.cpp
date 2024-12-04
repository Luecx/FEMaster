//
// Created by Finn Eggers on 05.09.23.
//
#include "../cos/cylindrical_system.h"
#include "../cos/rectangular_system.h"
#include "../model/beam/b23.h"
#include "../model/solid/c3d10.h"
#include "../model/solid/c3d15.h"
#include "../model/solid/c3d20.h"
#include "../model/solid/c3d20r.h"
#include "../model/solid/c3d4.h"
#include "../model/solid/c3d6.h"
#include "../model/solid/c3d8.h"
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
        } else {
            logging::warning(false, "Unknown command: ", m_current_line.line());
            next_line();
        }
    }
}

void Reader::process_nodes() {
    // read NODE
    // check if NSET is defined, if not, use NALL
    m_model->_data->node_sets.activate(m_current_line.parse<std::string>("NSET", "NALL"));

    while (next_line().type() == DATA_LINE) {
        int node_id = m_current_line.get_value(0, 0);
        Precision x = m_current_line.get_value(1, (Precision) 0);
        Precision y = m_current_line.get_value(2, (Precision) 0);
        Precision z = m_current_line.get_value(3, (Precision) 0);
        m_model->set_node(node_id, x, y, z);
    }
}

void Reader::process_elements() {
    // read ELEMENT
    // check if ELSET is defined, if not, use EALL
    m_model->_data->elem_sets.activate(m_current_line.parse<std::string>("ELSET", "EALL"));

    auto type = m_current_line.require<std::string>("TYPE");


    auto gather_values = [&](Index req_values) -> std::vector<ID> {
        std::vector<ID> values;
        for(Index i = 1; i < m_current_line.values().size(); i++){
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
        ID id = m_current_line.get_value(0, 0);
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
        } else if (type == "B23") {
            auto values = gather_values(2);
            m_model->set_element<fem::model::B23>(id, values[0], values[1]);
        }else {
            logging::warning(false, "Unknown element type ", type);
            return;
        }
    }
}

void Reader::process_surfaces() {
    // Activate the specified element set (if any) before processing surfaces.
    m_model->_data->surface_sets.activate(m_current_line.require<std::string>("SFSET", "NAME"));

    while(next_line().type() == DATA_LINE) {
        // Read the number of parts in the current line.
        int num_values = static_cast<int>(m_current_line.values().size());

        int id = -1; // Default value for ID if not specified.
        int elem_id = -1;
        int surf_side = -1;
        std::string set_name;

        if (num_values == 3) {
            // Three parts: ID, elementID, sideID
            id      = m_current_line.get_value(0, 0);
            elem_id = m_current_line.get_value(1, 0);
            if (m_current_line.values()[2][0] == 'S') {
                surf_side = std::stoi(m_current_line.values()[2].substr(1));
            } else {
                surf_side = m_current_line.get_value(2, 0);
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
                    surf_side = m_current_line.get_value(1, 0);
                }

                m_model->set_surface(id, elem_id, surf_side);
            } else {
                // First value is element set name
                set_name = first_value;
                if (m_current_line.values()[1][0] == 'S') {
                    surf_side = std::stoi(m_current_line.values()[1].substr(1));
                } else {
                    surf_side = m_current_line.get_value(1, 0);
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
    m_model->_data->node_sets.activate(setname);
    while (next_line().type() == DATA_LINE) {
        if (generate) {
            // require exactly 2 or 3 values, not more, not less
            logging::warning(m_current_line.count_values() == 2
                            || m_current_line.count_values() == 3, "GENERATE requires exactly 2 or 3 values.");

            ID id1 = m_current_line.get_value(0, 0);
            ID id2 = m_current_line.get_value(1, 0);
            ID inc = m_current_line.get_value(2, 1);
            for (ID i = id1; i <= id2; i += inc) {
                m_model->_data->node_sets.get()->add(i);
            }
        } else {
            for (const auto& id : m_current_line.values()) {
                m_model->_data->node_sets.get()->add(std::stoi(id));
            }
        }
    }
}

void Reader::process_elset() {
    // read ELSET, NAME=xyz
    bool generate = m_current_line.has_key("GENERATE");

    m_model->_data->elem_sets.activate(m_current_line.require<std::string>("NAME", "ELSET"));
    while (next_line().type() == DATA_LINE) {

        if (generate) {
            // require exactly 2 or 3 values, not more, not less
            logging::warning(m_current_line.count_values() == 2
                            || m_current_line.count_values() == 3, "GENERATE requires exactly 2 or 3 values.");

            ID id1 = m_current_line.get_value(0, 0);
            ID id2 = m_current_line.get_value(1, 0);
            ID inc = m_current_line.get_value(2, 1);
            for (ID i = id1; i <= id2; i += inc) {
                m_model->_data->elem_sets.get()->add(i);
            }
        } else {
            for (const auto& id : m_current_line.values()) {
                m_model->_data->elem_sets.get()->add(std::stoi(id));
            }
        }
    }
}

void Reader::process_sfset() {
    // read ELSET, NAME=xyz
    bool generate = m_current_line.has_key("GENERATE");
    m_model->_data->surface_sets.activate(m_current_line.require<std::string>("NAME", "SFSET"));

    while (next_line().type() == DATA_LINE) {
        if (generate) {
            // require exactly 2 or 3 values, not more, not less
            logging::warning(m_current_line.count_values() == 2
                            || m_current_line.count_values() == 3, "GENERATE requires exactly 2 or 3 values.");

            ID id1 = m_current_line.get_value(0, 0);
            ID id2 = m_current_line.get_value(1, 0);
            ID inc = m_current_line.get_value(2, 1);
            for (ID i = id1; i <= id2; i += inc) {
                m_model->_data->surface_sets.get()->add(i);
            }
        } else {
            for (const auto& id : m_current_line.values()) {
                m_model->_data->surface_sets.get()->add(std::stoi(id));
            }
        }
    }
}

void Reader::process_material() {
    // read MATERIAL, NAME=xyz
    m_model->_data->materials.activate(m_current_line.require<std::string>("NAME"));
    next_line();
}

void Reader::process_elastic() {
    // read ELASTIC, TYPE=xyz
    auto type = m_current_line.require<std::string>("TYPE");
    next_line();
    if (type == "ISO" || type == "ISOTROPIC") {
        auto E = m_current_line.get_value(0, 0.0f);
        auto n = m_current_line.get_value(1, 0.0f);
        m_model->_data->materials.get()->set_elasticity<fem::material::IsotropicElasticity>(E, n);
    }
    if (type == "ORTHO" || type == "ORTHOTROPIC") {
        auto E1 = m_current_line.get_value(0, 0.0f);
        auto E2 = m_current_line.get_value(1, 0.0f);
        auto E3 = m_current_line.get_value(2, 0.0f);
        auto G23 = m_current_line.get_value(3, 0.0f);
        auto G13 = m_current_line.get_value(4, 0.0f);
        auto G12 = m_current_line.get_value(5, 0.0f);
        auto nu23 = m_current_line.get_value(6, 0.0f);
        auto nu13 = m_current_line.get_value(7, 0.0f);
        auto nu12 = m_current_line.get_value(8, 0.0f);
        m_model->_data->materials.get()->set_elasticity<fem::material::OrthotropicElasticity>(
            E1, E2, E3, G23, G13, G12, nu23, nu13, nu12);
    }
    next_line();
}

void Reader::process_density() {
    // read DENSITY
    next_line();
    auto v = m_current_line.get_value(0, 0.0f);
    m_model->_data->materials.get()->set_density(v);
    next_line();
}

void Reader::process_thermal_expansion() {
    // read THERMAL EXPANSION
    next_line();
    auto v = m_current_line.get_value(0, 0.0f);
    m_model->_data->materials.get()->set_thermal_expansion(v);
    next_line();
}

void Reader::process_profile() {
    auto mat = m_current_line.require<std::string>("NAME", "PROFILE");
    next_line();
    Precision A = m_current_line.get_value(0, 0.0f);
    Precision Iy = m_current_line.get_value(1, 0.0f);
    Precision Iz = m_current_line.get_value(2, 0.0f);
    Precision Jt = m_current_line.get_value(3, 0.0f);
    m_model->_data->profiles.activate(mat, A, Iy, Iz, Jt);
    next_line();
}

void Reader::process_solid_section() {
    auto mat = m_current_line.require<std::string>("MAT", "MATERIAL");
    auto els = m_current_line.require<std::string>("ELSET");
    m_model->solid_section(els, mat);
    next_line();
}

void Reader::process_beam_section() {
    auto mat = m_current_line.require<std::string>("MAT", "MATERIAL");
    auto els = m_current_line.require<std::string>("ELSET");
    auto profile = m_current_line.require<std::string>("PROFILE");
    next_line();
    Precision n1x = m_current_line.get_value(0, 0.0f);
    Precision n1y = m_current_line.get_value(1, 0.0f);
    Precision n1z = m_current_line.get_value(2, 0.0f);
    Vec3 n1 = Vec3(n1x, n1y, n1z);
    m_model->beam_section(els, mat, profile, n1);
    next_line();
}

void Reader::process_cload() {
    // read CLOAD, LOAD_COLLECTOR=xyz
    // NSET, lx, ly, lz
    // id, lx, ly, lz
    // ...
    m_model->_data->load_cols.activate(m_current_line.require<std::string>("LOAD_COLLECTOR"));
    while (next_line().type() == DATA_LINE) {
        auto str = m_current_line.values()[0];
        auto lx  = m_current_line.get_value(1, 0.0f);
        auto ly  = m_current_line.get_value(2, 0.0f);
        auto lz  = m_current_line.get_value(3, 0.0f);

        auto mx  = m_current_line.get_value(4, 0.0f);
        auto my  = m_current_line.get_value(5, 0.0f);
        auto mz  = m_current_line.get_value(6, 0.0f);

        if (m_model->_data->node_sets.has(str)) {
            m_model->add_cload(str, Vec6(lx, ly, lz, mx, my, mz));
        } else {
            m_model->add_cload(std::stoi(m_current_line.values()[0]), Vec6(lx, ly, lz, mx, my, mz));
        }
    }
}

void Reader::process_dload() {
    // read DLOAD, LOAD_COLLECTOR=xyz
    // SFSET, lx, ly, lz
    // id, lx, ly, lz
    // ...
    m_model->_data->load_cols.activate(m_current_line.require<std::string>("LOAD_COLLECTOR"));

    while (next_line().type() == DATA_LINE) {
        auto str = m_current_line.values()[0];
        auto lx  = m_current_line.get_value(1, 0.0f);
        auto ly  = m_current_line.get_value(2, 0.0f);
        auto lz  = m_current_line.values().size() > 3 ? std::stof(m_current_line.values()[3]) : 0;

        if (m_model->_data->surface_sets.has(str)) {
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
    m_model->_data->load_cols.activate(m_current_line.require<std::string>("LOAD_COLLECTOR"));
    while (next_line().type() == DATA_LINE) {
        auto str = m_current_line.values()[0];
        auto lx  = m_current_line.get_value(1, 0.0f);
        auto ly  = m_current_line.get_value(2, 0.0f);
        auto lz  = m_current_line.get_value(3, 0.0f);

        if (m_model->_data->elem_sets.has(str)) {
            m_model->add_vload(str, Vec3(lx, ly, lz));
        } else {
            m_model->add_vload(std::stoi(m_current_line.values()[0]), Vec3(lx, ly, lz));
        }
    }
}

void Reader::process_tload() {
    // read TLOAD, LOAD_COLLECTOR=xyz, temperature field=xyz, reference temperature=xyz
    auto lod_col  = m_current_line.require<std::string>("LOAD_COLLECTOR");
    auto temp_fi  = m_current_line.require<std::string>("TEMPERATUREFIELD");
    auto ref_temp = m_current_line.require<Precision  >("REFERENCETEMPERATURE");

    m_model->_data->load_cols.activate(lod_col);
    m_model->add_tload(temp_fi, ref_temp);
    next_line();
}

void Reader::process_support() {
    // read SUPPORT, SUPPORT_COLLECTOR=xyz
    // NSET, lx, ly, lz
    // id, lx, ly, lz
    // ...
    // use empty field if no support
    auto supp_col    = m_current_line.require<std::string>("SUPPORT_COLLECTOR");
    auto orientation = m_current_line.parse<std::string>("ORIENTATION", "");

    m_model->_data->supp_cols.activate(supp_col);

    while (next_line().type() == DATA_LINE) {
        auto str = m_current_line.values()[0];
        StaticVector<6> constraint;

        for(Dim i = 0; i < 6; i++) {
            bool is_given = i < m_current_line.values().size() - 1;
            bool is_empty = !is_given || m_current_line.values()[i + 1].empty();

            if (is_given && !is_empty) {
                constraint(i) = std::stof(m_current_line.values()[i + 1]);
            } else {
                constraint(i) = NAN;
            }
        }

        if (m_model->_data->node_sets.has(str)) {
            m_model->add_support(str, constraint, orientation);
        } else {
            m_model->add_support(std::stoi(str), constraint, orientation);
        }
    }
}

void Reader::process_temperature() {
    // read TEMPERATURE, NAME=xyz
    // NSET, value
    // id, value
    // ...
    auto name = m_current_line.require<std::string>("NAME");

    while (next_line().type() == DATA_LINE) {
        auto str   = m_current_line.values()[0];
        auto value = m_current_line.get_value(1, 0.0f);

        if (m_model->_data->node_sets.has(str)) {
            for (auto id : *m_model->_data->node_sets.get(str)) {
                m_model->set_field_temperature(name, id, value);
            }
        } else {
            m_model->set_field_temperature(name, std::stoi(str), value);
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
            auto x = m_current_line.get_value(0, 0.0f);
            auto y = m_current_line.get_value(1, 0.0f);
            auto z = m_current_line.get_value(2, 0.0f);
            m_model->add_coordinate_system<cos::RectangularSystem>(name, Vec3{x,y,z});
        } else if (m_current_line.count_values() == 6) {
            auto x1 = m_current_line.get_value(0, 0.0f);
            auto y1 = m_current_line.get_value(1, 0.0f);
            auto z1 = m_current_line.get_value(2, 0.0f);
            auto x2 = m_current_line.get_value(3, 0.0f);
            auto y2 = m_current_line.get_value(4, 0.0f);
            auto z2 = m_current_line.get_value(5, 0.0f);
            m_model->add_coordinate_system<cos::RectangularSystem>(name, Vec3{x1,y1,z1}, Vec3{x2,y2,z2});
        } else if (m_current_line.count_values() == 9) {
            auto x1 = m_current_line.get_value(0, 0.0f);
            auto y1 = m_current_line.get_value(1, 0.0f);
            auto z1 = m_current_line.get_value(2, 0.0f);
            auto x2 = m_current_line.get_value(3, 0.0f);
            auto y2 = m_current_line.get_value(4, 0.0f);
            auto z2 = m_current_line.get_value(5, 0.0f);
            auto x3 = m_current_line.get_value(6, 0.0f);
            auto y3 = m_current_line.get_value(7, 0.0f);
            auto z3 = m_current_line.get_value(8, 0.0f);
            m_model->add_coordinate_system<cos::RectangularSystem>(name, Vec3{x1,y1,z1}, Vec3{x2,y2,z2}, Vec3{x3,y3,z3});
        } else {
            logging::error(false, "Cannot create coordinate system with", m_current_line.count_values(),
                " values. Needs to be 3, 6 or 9.");
        }
    } else if (type== "CYLINDRICAL") {
        if (m_current_line.count_values() == 9) {
            auto x1 = m_current_line.get_value(0, 0.0f);
            auto y1 = m_current_line.get_value(1, 0.0f);
            auto z1 = m_current_line.get_value(2, 0.0f);
            auto x2 = m_current_line.get_value(3, 0.0f);
            auto y2 = m_current_line.get_value(4, 0.0f);
            auto z2 = m_current_line.get_value(5, 0.0f);
            auto x3 = m_current_line.get_value(6, 0.0f);
            auto y3 = m_current_line.get_value(7, 0.0f);
            auto z3 = m_current_line.get_value(8, 0.0f);
            m_model->add_coordinate_system<cos::CylindricalSystem>(name, Vec3{x1,y1,z1}, Vec3{x2,y2,z2}, Vec3{x3,y3,z3});
        } else {
            logging::error(false, "Cannot create coordinate system with", m_current_line.count_values(),
                " values. Needs to be exactly 9.");
        }
    }

    else {
        logging::error(false, "Unknown coordinate system type: ", type);
    }
    next_line();
}

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
    } else {
        logging::error(false, "Unknown connector type: ", type);
    }

    m_model->add_connector(nset1, nset2, coord, ctype);

    next_line();
}


void Reader::process_coupling() {
    // read COUPLING, MASTER=xyz, SLAVE=xyz, DOFS=xyz, TYPE=xyz
    auto master_set = m_current_line.require<std::string>("MASTER");
    auto type       = m_current_line.require<std::string>("TYPE");

    // either SLAVE or SURFACE can be defined. slave refers to a node set and surface to a surface set
    auto slave_set  = m_current_line.parse  <std::string>("SLAVE", "");
    auto surface    = m_current_line.parse  <std::string>("SFSET", "");

    next_line();

    Dofs dof_mask {false, false, false, false, false, false};
    for (Dim i = 0; i < 6; i++) {
        if (i < m_current_line.values().size()) {
            dof_mask(i) = std::stof(m_current_line.values()[i]) > 0;
        }
    }

    if (type == "KINEMATIC") {
        if (surface.empty()) {
            m_model->add_coupling(master_set, slave_set, dof_mask, constraint::CouplingType::KINEMATIC, false);
        } else {
            m_model->add_coupling(master_set, surface, dof_mask, constraint::CouplingType::KINEMATIC, true);
        }
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
        } else if (m_current_line.command() == "ORIENTATION") {
            process_loadcase_linear_static_topo_orientation(&lc);
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
        auto id         = m_current_line.get_value(0, 0);
        auto ds         = m_current_line.get_value(1, 0.0f);

        lc->density(id) = ds;
    }
}
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
void Reader::process_loadcase_linear_static_topo_exponent(fem::loadcase::LinearStaticTopo* lc) {
    next_line();
    auto exp     = m_current_line.get_value(0, 0.0f);
    lc->exponent = exp;
    next_line();
}

}    // namespace fem::reader