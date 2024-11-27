//
// Created by Finn Eggers on 03.09.23.
//

#include "model.h"
#include "../mattools/numerate_dofs.h"
#include "../mattools/assemble.h"
#include "../mattools/assemble_bc.h"
#include "../mattools/lump_matrix.h"

namespace fem {

namespace model {

void Model::assign_sections() {
    for (Section::Ptr ptr: this->_data->sections) {
        for (ID elem_id = 0; elem_id < this->_data->max_elems; elem_id++) {
            if (this->_data->elements[elem_id] != nullptr) {
                this->_data->elements[elem_id]->set_section(ptr);
            }
        }
    }
}


SystemDofIds Model::build_unconstrained_index_matrix() {
    // first build a boolean matrix of which dofs are present in the system
    SystemDofs mask{this->_data->max_nodes, 6};
    mask.fill(false);

    // go through each elements and mask the dofs for all the nodes
    for (auto &e: _data->elements) {
        if (e != nullptr) {
            auto dofs = e->dofs();
            for (ID node_local_id = 0; node_local_id < e->n_nodes(); node_local_id++) {
                ID node_id = e->nodes()[node_local_id];
                for (ID dof = 0; dof < 6; dof++) {
                    mask(node_id, dof) |= dofs(0, dof);
                }
            }
        }
    }

    // go through all _couplings and mask the master dof
    for (auto &c: this->_data->couplings) {
        ID master_id = c.master_node;
        auto master_dofs = c.master_dofs(mask, *_data);
        for (ID dof = 0; dof < 6; dof++) {
            mask(master_id, dof) |= master_dofs(0, dof);
        }
    }

    auto res = mattools::numerate_dofs(mask);

    return res;
}

NodeData Model::build_support_matrix(std::vector<std::string> supp_sets) {
    NodeData disp_matrix{this->_data->max_nodes, 6};
    disp_matrix.fill(std::numeric_limits<Precision>::quiet_NaN());

    for (auto &key: supp_sets) {
        auto data = this->_support_sets.get(key);
        mattools::assemble_bc(disp_matrix, *data, mattools::DuplicateHandling::SET);
    }
    return disp_matrix;
}

NodeData Model::build_load_matrix(std::vector<std::string> load_sets) {
    NodeData load_matrix{this->_data->max_nodes, 6};
    load_matrix.setZero();

    for (auto &key: load_sets) {
        auto data = this->_load_sets.get(key);
        mattools::assemble_bc(load_matrix, *data, mattools::DuplicateHandling::ADD);
    }
    return load_matrix;
}

SparseMatrix Model::build_constraint_matrix   (SystemDofIds& indices, Precision characteristic_stiffness) {
    constraint::Equations equations;

    for (auto &c: this->_data->couplings) {
        auto coupling_equations = c.get_equations(indices, *_data);
        for (auto &eq: coupling_equations) {
            equations.push_back(eq);
        }
    }

    for (auto &t: this->_data->ties) {
        auto tie_equations = t.get_equations(indices, *_data);
        for (auto &eq: tie_equations) {
            equations.push_back(eq);
        }
    }

    for (auto &t : this->_data->connectors) {
        auto connector_equations = t.get_equations(indices, *_data);
        for (auto &eq: connector_equations) {
            equations.push_back(eq);
        }
    }

    for (auto &eq: this->_data->equations) {
        equations.push_back(eq);
    }

    if (equations.empty()) {
        return SparseMatrix{0, indices.maxCoeff() + 1};
    }
    auto triplets = constraint::Equation::get_triplets(equations, indices, 0);

    SparseMatrix matrix{(Eigen::Index) equations.size(), indices.maxCoeff() + 1};
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix *= characteristic_stiffness;
    return matrix;
}



SparseMatrix Model::build_stiffness_matrix(SystemDofIds &indices, ElementData stiffness_scalar) {
    auto lambda = [&](const ElementPtr &el, Precision* storage) {
        if (auto sel = el->as<StructuralElement>()) {
            MapMatrix stiff = sel->stiffness(storage);
            stiff *= (stiffness_scalar.rows() > 0) ? stiffness_scalar(sel->elem_id, 0) : 1.0;
            return stiff;
        } else {
            MapMatrix mat{storage, 0, 0};
            return mat;
        }
    };
    auto res = mattools::assemble_matrix(_data->elements, indices, lambda);
    return res;
}

SparseMatrix Model::build_lumped_mass_matrix(SystemDofIds& indices) {
    auto lambda = [&](const ElementPtr &el, Precision* storage) {
        if (auto sel = el->as<StructuralElement>()) {
            MapMatrix element_mass = sel->mass(storage);
            return element_mass;
        } else {
            MapMatrix mat{storage, 0, 0};
            return mat;
        }
    };
    auto res = mattools::assemble_matrix(_data->elements, indices, lambda);
    return res;
}


}
}