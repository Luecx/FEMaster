#include <memory>

#include "line.h"
#include "file.h"
#include "../model/model.h"
#include "../loadcase/linear_static.h"
#include "../loadcase/linear_static_topo.h"

namespace fem::reader {

class Reader {
    private:
    struct ProblemData {
        int highest_node_id = -1;
        int highest_element_id = -1;
        int current_loadcase_num = 1;
    };

    File        m_file;
    Writer      m_writer;
    Line        m_current_line;
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


    void process();
    void process_nodes();
    void process_elements();
    void process_nset();
    void process_elset() ;
    void process_material() ;
    void process_elastic() ;
    void process_density();
    void process_solid_section();
    void process_cload();
    void process_support();

    void process_loadcase();
    void process_loadcase_linear_static();
    void process_loadcase_linear_static_topo();
    void process_loadcase_linear_static_support (fem::loadcase::LinearStatic* lc);
    void process_loadcase_linear_static_load    (fem::loadcase::LinearStatic* lc);
    void process_loadcase_linear_static_solver  (fem::loadcase::LinearStatic* lc);

    void process_loadcase_linear_static_topo_density (fem::loadcase::LinearStaticTopo* lc);
    void process_loadcase_linear_static_topo_exponent(fem::loadcase::LinearStaticTopo* lc);



    public:
    Reader(const std::string& file_path) :
            m_file(file_path),
            m_file_path(file_path),
            m_writer(file_path + ".res") {}

    void read() {
        m_file   = File(m_file_path);
        m_writer = Writer(m_file_path + ".res");

        // First stage
        analyse();

        // create the model
        m_model = std::make_shared<fem::model::Model>(m_data.highest_node_id + 1,
                                                      m_data.highest_element_id + 1);

        // Reset File for the next stage
        m_file = File(m_file_path);

        // Second stage
        process();

        // close the writer
        m_writer.close();
    }
};

}
