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

SystemDofIds Model::build_unconstrained_index_matrix() {
    // first build a boolean matrix of which dofs are present in the system
    SystemDofs mask{this->max_nodes, 6};
    mask.fill(false);

    // go through each elements and mask the dofs for all the nodes
    for (auto &e: elements) {
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

    // go through all couplings and mask the master dof
    for (auto &c: this->couplings) {
        ID master_id = c.master_node;
        auto master_dofs = c.master_dofs(mask);
        for (ID dof = 0; dof < 6; dof++) {
            mask(master_id, dof) |= master_dofs(0, dof);
        }
    }

    auto res = mattools::numerate_dofs(mask);

    return res;
}

NodeData Model::build_support_matrix(std::vector<std::string> supp_sets) {
    NodeData disp_matrix{this->max_nodes, 6};
    disp_matrix.fill(std::numeric_limits<Precision>::quiet_NaN());

    for (auto &key: supp_sets) {
        NodeData &data = this->support_sets.get(key);
        mattools::assemble_bc(disp_matrix, data, mattools::DuplicateHandling::SET);
    }
    return disp_matrix;
}

NodeData Model::build_load_matrix(std::vector<std::string> load_sets) {
    NodeData load_matrix{this->max_nodes, 6};
    load_matrix.setZero();

    for (auto &key: load_sets) {
        NodeData &data = this->load_sets.get(key);
        mattools::assemble_bc(load_matrix, data, mattools::DuplicateHandling::ADD);
    }
    return load_matrix;
}

SparseMatrix Model::build_constraint_matrix   (SystemDofIds& indices) {
    TripletList triplets;
    int rows = 0;
    for (auto &c: this->couplings) {
        auto coupling_triplets = c.get_equations(indices, node_coords, rows);
        for (auto &t: coupling_triplets) {
            triplets.push_back(t);
            rows = std::max(rows, t.row());
        }
    }

    if (rows == 0) {
        return SparseMatrix{0, indices.maxCoeff() + 1};
    }

    SparseMatrix matrix{rows + 1, indices.maxCoeff() + 1};
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return matrix;
}



SparseMatrix Model::build_stiffness_matrix(SystemDofIds &indices, ElementData stiffness_scalar) {
    auto lambda = [&](const ElementPtr &el, Precision* storage) {
        MapMatrix stiff = el->stiffness(node_coords, storage);
        stiff *= (stiffness_scalar.rows() > 0) ? stiffness_scalar(el->elem_id, 0) : 1.0;
        return stiff;
    };
    auto res = mattools::assemble_matrix(elements, indices, lambda);
    return res;
}

SparseMatrix Model::build_lumped_mass_matrix(SystemDofIds& indices) {
    auto lambda = [&](const ElementPtr &el, Precision* storage) {
        MapMatrix element_mass = el->mass(node_coords, storage);
        return element_mass;
    };
    auto res = mattools::assemble_matrix(elements, indices, lambda);
    return res;
}


}
}