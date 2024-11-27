
#include "support.h"

#include "../model/element/element.h"
#include "../model/model_data.h"

namespace fem {

using namespace model;

void apply_to_node(model::ModelData& model_data, NodeData& bc, constraint::Equations& equations,
    ID node_id, Vec6& values, const cos::CoordinateSystem::Ptr& coordinate_system) {

    Vec3 position = model_data.get(POSITION).row(node_id);
    for(int i = 0; i < 6; i++) {
        if (!std::isnan(values[i])) {

            // check if transformation is needed, if so, add equations inside of adding to the bc
            if (coordinate_system != nullptr) {

                // get the transformation matrix
                auto transformation = coordinate_system->get_axes(position);

                // values
                Vec3 vals = transformation.col(i % 6);

                constraint::EquationEntry en1 = {node_id, (Dim) ((i / 3) * 3), vals[0]};
                constraint::EquationEntry en2 = {node_id, (Dim) ((i / 3) * 3 + 1), vals[1]};
                constraint::EquationEntry en3 = {node_id, (Dim) ((i / 3) * 3 + 2), vals[2]};
                equations.push_back(constraint::Equation({en1, en2, en3}));
            }

            else {
                bc(node_id, i) = values[i];
            }
        }
    }
}

void NSupport::apply(model::ModelData& model_data, NodeData& bc, constraint::Equations& equations){
    for(ID id: *(this->region)) {
        apply_to_node(model_data, bc, equations, id, this->values, this->coordinate_system);
    }
}

void ESupport::apply(model::ModelData& model_data, NodeData& bc, constraint::Equations& equations){
    for(ID el_id: *(this->region)) {
        if (model_data.elements[el_id] == nullptr)
            continue;

        auto el = model_data.elements[el_id];
        for(int local_idx = 0; local_idx < el->n_nodes(); local_idx++) {
            apply_to_node(model_data, bc, equations, el->nodes()[local_idx], this->values, this->coordinate_system);
        }
    }
}

void SSupport::apply(model::ModelData& model_data, NodeData& bc, constraint::Equations& equations){
    for(ID el_id: *(this->region)) {
        if (model_data.surfaces[el_id] == nullptr)
            continue;

        auto el = model_data.surfaces[el_id];
        for(Index local_idx = 0; local_idx < el->n_nodes; local_idx++) {
            apply_to_node(model_data, bc, equations, el->nodes()[local_idx], this->values, this->coordinate_system);
        }
    }
}

}