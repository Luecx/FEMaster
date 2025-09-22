#include "../model/beam/b33.h"
#include "../model/truss/truss.h"
#include "../model/model.h"
#include "../model/pointelem/point.h"
#include "../model/shell/s3.h"
#include "../model/shell/s4.h"
#include "../model/shell/s4_mitc.h"
#include "../model/shell/s6.h"
#include "../model/shell/s8.h"
#include "../model/solid/c3d10.h"
#include "../model/solid/c3d15.h"
#include "../model/solid/c3d20.h"
#include "../model/solid/c3d20r.h"
#include "../model/solid/c3d4.h"
#include "../model/solid/c3d6.h"
#include "../model/solid/c3d8.h"
#include "reader.h"

namespace fem::reader {
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
        } else if (type == "B33") {
            auto values = gather_values(2);
            m_model->set_element<fem::model::B33>(id, values[0], values[1]);
        } else if (type == "T3") {
            auto values = gather_values(2);
            m_model->set_element<fem::model::T3>(id, values[0], values[1]);
        } else if (type == "P") {
            auto values = gather_values(1);
            m_model->set_element<fem::model::Point>(id, values[0]);
        } else if (type == "S3") {
            auto values = gather_values(3);
            m_model->set_element<fem::model::S3>(id, values[0], values[1], values[2]);
        } else if (type == "S4") {
            auto values = gather_values(4);
            m_model->set_element<fem::model::S4>(id, values[0], values[1], values[2], values[3]);
        } else if (type == "MITC4") {
            logging::error(false, "Element type MITC4 is supported but shall not be used yet");
            auto values = gather_values(4);
            m_model->set_element<fem::model::MITC4>(id, values[0], values[1], values[2], values[3]);
        } else if (type == "S6") {
            auto values = gather_values(6);
            m_model->set_element<fem::model::S6>(id, values[0], values[1], values[2], values[3], values[4], values[5]);
        } else if (type == "S8") {
            auto values = gather_values(8);
            m_model->set_element<fem::model::S8>(id, values[0], values[1], values[2], values[3],
                                                     values[4], values[5], values[6], values[7]);
        } else {
            logging::warning(false, "Unknown element type ", type);
            return;
        }
    }
}
} // namespace fem::reader
