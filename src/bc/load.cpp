//
// Created by Luecx on 28.11.2024.
//

#include "../model/model_data.h"

#include "../model/element/element_structural.h"
#include "../model/element/element.h"

#include "load.h"


void fem::CLoad::apply(fem::model::ModelData& model_data, fem::NodeData& bc) {
    (void) model_data;
    for (auto& node_id : *region) {
        for (Dim i = 0; i < 6; i++) {
            if (!std::isnan(values[i])) {
                bc(node_id, i) += values[i];
            }
        }
    }
}

void fem::DLoad::apply(fem::model::ModelData& model_data, fem::NodeData& bc) {
    auto& node_positions = model_data.get(model::POSITION);
    for(auto& surf_id : *region) {
        model_data.surfaces[surf_id]->apply_dload(node_positions, bc, values);
    }
}

void fem::PLoad::apply(fem::model::ModelData& model_data, fem::NodeData& bc) {
    auto& node_positions = model_data.get(model::POSITION);
    for(auto& surf_id : *region) {
        model_data.surfaces[surf_id]->apply_pload(node_positions, bc, pressure);
    }
}

void fem::VLoad::apply(fem::model::ModelData& model_data, fem::NodeData& bc) {
    for(auto& el_id : *region) {
        auto& el_ptr = model_data.elements[el_id];
        auto els_ptr = el_ptr->as<model::StructuralElement>();
        if (els_ptr) {
            els_ptr->apply_vload(bc, values);
        }
    }
}

void fem::TLoad::apply(fem::model::ModelData& model_data, fem::NodeData& bc) {
    for (auto& el_ptr: model_data.elements) {
        if (auto sel = el_ptr->as<model::StructuralElement>()) {
            sel->apply_tload(bc, *temp_field, ref_temp);
        }
    }
}
