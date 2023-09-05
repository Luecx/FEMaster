#include <memory>

#include "line.h"
#include "file.h"
#include "../model/model.h"

namespace fem::reader {

class Reader {
    private:
    struct ProblemData {
        int highest_node_id = -1;
        int highest_element_id = -1;
    };

    File m_file;
    Line m_current_line;
    ProblemData m_data;
    std::string m_file_path;

    std::shared_ptr<fem::model::Model> m_model;

    Line& next_line() {
        m_current_line = m_file.next();
        return m_current_line;
    }

    void analyse() {
        // read first line
        next_line();

        while(m_current_line.type() != END_OF_FILE){
            while(m_current_line.type() != KEYWORD_LINE && m_current_line.type() != END_OF_FILE){
                next_line();
            }

            if (m_current_line.command() == "NODE") {
                analyse_nodes();
            } else if (m_current_line.command() == "ELEMENT") {
                analyse_elements();
            } else{
                next_line();
            }
        }
    }

    void analyse_nodes() {
        while (next_line().type() == DATA_LINE) {
            int node_id = std::stoi(m_current_line.values()[0]);
            m_data.highest_node_id = std::max(m_data.highest_node_id, node_id);
        }
    }

    void analyse_elements() {
        while (next_line().type() == DATA_LINE) {
            int element_id = std::stoi(m_current_line.values()[0]);
            m_data.highest_element_id = std::max(m_data.highest_element_id, element_id);
        }
    }


    void process() {
        // read first line
        next_line();

        while(m_current_line.type() != END_OF_FILE){
            while(m_current_line.type() != KEYWORD_LINE && m_current_line.type() != END_OF_FILE){
                next_line();
            }
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
            }
        }
    }

    void process_nodes() {
        // read NODE
        // check if NSET is defined, if not, use NALL
        m_model->activate_node_set(m_current_line.parse<std::string>("NSET", ""));

        while (next_line().type() == DATA_LINE) {
            int node_id = std::stoi(m_current_line.values()[0]);
            int x       = std::stof(m_current_line.values()[1]);
            int y       = std::stof(m_current_line.values()[2]);
            int z       = m_current_line.values().size() > 3 ? std::stof(m_current_line.values()[3]) : 0;
            m_model->set_node(node_id, x, y, z);
        }
    }

    void process_elements() {
        // read ELEMENT
        // check if ELSET is defined, if not, use EALL
        m_model->activate_element_set(m_current_line.parse<std::string>("ELSET", ""));

        auto type = m_current_line.require<std::string>("TYPE");

        while (next_line().type() == DATA_LINE) {
            int id = std::stoi(m_current_line.values()[0]);
            if (type == "C3D8"){
                m_model->set_element<fem::model::C3D8> (id
                                         , std::stoi(m_current_line.values()[1])
                                         , std::stoi(m_current_line.values()[2])
                                         , std::stoi(m_current_line.values()[3])
                                         , std::stoi(m_current_line.values()[4])
                                         , std::stoi(m_current_line.values()[5])
                                         , std::stoi(m_current_line.values()[6])
                                         , std::stoi(m_current_line.values()[7])
                                         , std::stoi(m_current_line.values()[8]));
            }
            if (type == "C3D20"){
                m_model->set_element<fem::model::C3D20> (id
                                         , std::stoi(m_current_line.values()[1])
                                         , std::stoi(m_current_line.values()[2])
                                         , std::stoi(m_current_line.values()[3])
                                         , std::stoi(m_current_line.values()[4])
                                         , std::stoi(m_current_line.values()[5])
                                         , std::stoi(m_current_line.values()[6])
                                         , std::stoi(m_current_line.values()[7])
                                         , std::stoi(m_current_line.values()[8])
                                         , std::stoi(m_current_line.values()[9])
                                         , std::stoi(m_current_line.values()[10])
                                         , std::stoi(m_current_line.values()[11])
                                         , std::stoi(m_current_line.values()[12])
                                         , std::stoi(m_current_line.values()[13])
                                         , std::stoi(m_current_line.values()[14])
                                         , std::stoi(m_current_line.values()[15])
                                         , std::stoi(m_current_line.values()[16])
                                         , std::stoi(m_current_line.values()[17])
                                         , std::stoi(m_current_line.values()[18]));
            }
        }
    }

    void process_nset() {
        // read NSET, NAME=xyz
        m_model->activate_node_set(m_current_line.require<std::string>("NAME"));
        while (next_line().type() == DATA_LINE) {
            for(const auto& id : m_current_line.values()) {
                m_model->active_nodeset().push_back(std::stoi(id));
            }
        }
    }

    void process_elset() {
        // read ELSET, NAME=xyz
        m_model->activate_element_set(m_current_line.require<std::string>("NAME"));
        while (next_line().type() == DATA_LINE) {
            for(const auto& id : m_current_line.values()) {
                m_model->active_elemset().push_back(std::stoi(id));
            }
        }
    }

    void process_material() {
        // read MATERIAL, NAME=xyz
        m_model->activate_material(m_current_line.require<std::string>("NAME"));
        next_line();
    }

    void process_elastic() {
        // read ELASTIC, TYPE=xyz
        auto type = m_current_line.require<std::string>("TYPE");
        next_line();
        if (type == "ISO"){
            auto E = std::stof(m_current_line.values()[0]);
            auto n = std::stof(m_current_line.values()[1]);
            m_model->active_material().set_elasticity<fem::material::IsotropicElasticity>(E, n);
        }
        next_line();
    }

    void process_density() {
        // read DENSITY
        next_line();
        auto v = std::stof(m_current_line.values()[0]);
        m_model->active_material().set_density(v);
        next_line();
    }

    public:
    Reader(const std::string& file_path) : m_file(file_path), m_file_path(file_path) {}

    void read() {
        m_file = File(m_file_path);

        // First stage
        analyse();

        // create the model
        m_model = std::make_shared<fem::model::Model>(m_data.highest_node_id + 1,
                                                      m_data.highest_element_id + 1);

        // Reset File for the next stage
        m_file = File(m_file_path);

        // Second stage
        process();

        std::cout<< *m_model << std::endl;
    }
};

}
