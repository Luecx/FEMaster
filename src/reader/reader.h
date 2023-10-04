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

    Line& next_line();

    void analyse();
    void analyse_nodes();
    void analyse_elements();

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
    void process_vload();
    void process_support();

    void process_loadcase();
    void process_loadcase_linear_static();
    void process_loadcase_linear_static_topo();
    void process_loadcase_linear_static_support     (fem::loadcase::LinearStatic* lc);
    void process_loadcase_linear_static_load        (fem::loadcase::LinearStatic* lc);
    void process_loadcase_linear_static_solver      (fem::loadcase::LinearStatic* lc);

    void process_loadcase_linear_static_topo_density (fem::loadcase::LinearStaticTopo* lc);
    void process_loadcase_linear_static_topo_exponent(fem::loadcase::LinearStaticTopo* lc);

    public:
    Reader(const std::string& file_path);

    void read();
};

}
