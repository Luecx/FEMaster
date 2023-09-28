//
// Created by Finn Eggers on 03.09.23.
//

#include "model.h"

namespace fem {

namespace model {

IndexMatrix Model::build_unconstrained_index_matrix() {

    // first build a boolean matrix of which dofs are present in the system
    BooleanMatrix mask{this->max_nodes, 6};
    mask.fill(false);

    // go through each elements and mask the dofs for all the nodes
    for (auto &e: elements) {
        if (e != nullptr) {
            auto dofs = e->dofs();
            for (ID node_local_id = 0; node_local_id < e->n_nodes(); node_local_id++) {
                ID node_id = e->nodes()[node_local_id];
                for (ID dof = 0; dof < 6; dof++) {
                    mask(node_id, dof) |= dofs(dof);
                }
            }
        }
    }

    // numerate the entries in the mask
    IndexMatrix res{this->max_nodes, 6};

    ID c = 0;
    for (ID n = 0; n < this->max_nodes; n++) {
        for (ID d = 0; d < 6; d++) {
            if (mask(n, d)) {
                res(n, d) = c++;
            } else {
                res(n, d) = -1;
            }
        }
    }

    return res;
}

IndexMatrix Model::build_constrained_index_matrix(DynamicVector &support_vector) {

    // first build a boolean matrix of which dofs are present in the system
    BooleanMatrix mask{this->max_nodes, 6};
    mask.fill(false);

    // go through each element and mask the dofs for all the nodes
    for (auto &e: elements) {
        if (e != nullptr) {
            auto dofs = e->dofs();
            for (ID node_local_id = 0; node_local_id < e->n_nodes(); node_local_id++) {
                ID node_id = e->nodes()[node_local_id];
                for (ID dof = 0; dof < 6; dof++) {
                    mask(node_id, dof) |= dofs(dof);
                }
            }
        }
    }

    // numerate the entries in the mask
    IndexMatrix res{this->max_nodes, 6};

    ID c = 0;
    ID c_unconstrained = 0;
    for (ID n = 0; n < this->max_nodes; n++) {
        for (ID d = 0; d < 6; d++) {
            if (mask(n, d)) {
                if (std::isnan(support_vector(c_unconstrained))) {
                    res(n, d) = c++;
                } else {
                    res(n, d) = -1;
                }
                c_unconstrained++;
            } else {
                res(n, d) = -1;
            }
        }
    }

    return res;
}


DynamicVector Model::build_support_vector(IndexMatrix &indices, std::vector <std::string> supp_sets) {
    DynamicVector support_vector{indices.maxCoeff() + 1};
    support_vector.fill(std::numeric_limits<Precision>::quiet_NaN());

    for (auto key: supp_sets) {
        NodeData &data = this->support_sets.get(key);

        for (int m = 0; m < data.rows(); m++) {
            for (int n = 0; n < data.cols(); n++) {
                auto idx = indices(m, n);
                if (idx >= 0 && !std::isnan(data(m, n))) {
                    logging::error(std::isnan(support_vector(idx)) || support_vector(idx) == data(m, n),
                                   "two different support collectors seem to both constrain the same "
                                   "node with different values, this is not allowed");
                    support_vector(idx) = data(m, n);
                }
            }
        }
    }
    return support_vector;
}

DynamicVector Model::build_load_vector(IndexMatrix &indices, std::vector <std::string> load_sets) {
    DynamicVector load_vector{indices.maxCoeff() + 1};
    load_vector.setZero();
    for (auto key: load_sets) {
        NodeData &data = this->load_sets.get(key);

        for (int m = 0; m < data.rows(); m++) {
            for (int n = 0; n < data.cols(); n++) {
                auto idx = indices(m, n);
                if (idx >= 0) {
                    load_vector(idx) += data(m, n);
                }
            }
        }
    }
    return load_vector;
}

SparseMatrix Model::build_stiffness_matrix(IndexMatrix &indices, ElementData stiffness_scalar) {
    constexpr int BATCH_SIZE = 1024 * 1024 * 16;
    constexpr int MAX_LOCAL_SIZE = 128;

    SparseMatrix matrix{indices.maxCoeff() + 1, indices.maxCoeff() + 1};
    SparseMatrixBuilder tripplets{};
    tripplets.reserve(BATCH_SIZE);

    alignas(64)
    Precision local_storage[MAX_LOCAL_SIZE * MAX_LOCAL_SIZE]{};

    for (auto el: elements) {
        if (el == nullptr) continue;

        // get the scalar if one defined
        Precision scalar = 1;
        if(stiffness_scalar.rows() > el->elem_id){
            scalar = stiffness_scalar(el->elem_id, 0);
        }

        // get the stiffness matrix and do remaining stuff
        MapMatrix K = el->stiffness(node_coords, local_storage);

        auto el_nodes = el->n_nodes();     // nodes in this element
        auto el_ndofs = K.rows() / el_nodes; // dofs per node

        auto el_dof_m = el->dofs();
        IndexVector dof_index_lookup{6}; // lookup to know that dof 2 for example equals index 5 or so

        int current_idx = 0;
        for (int i = 0; i < 6; ++i) {
            if (el_dof_m[i]) {
                dof_index_lookup[current_idx] = i;
                current_idx++;
            }
        }

        for (int i = 0; i < el_nodes; i++) {
            for (int j = 0; j < el_nodes; j++) {
                for (int idof = 0; idof < el_ndofs; idof++) {
                    for (int jdof = 0; jdof < el_ndofs; jdof++) {
                        int idof_idx = dof_index_lookup(idof);
                        int jdof_idx = dof_index_lookup(jdof);
                        ID global_i = indices(el->nodes()[i], idof_idx);
                        ID global_j = indices(el->nodes()[j], jdof_idx);

                        tripplets.emplace_back(global_i, global_j, scalar * K(i * el_ndofs + idof, j * el_ndofs + jdof));
                    }
                }
            }
        }

        logging::info(el->elem_id % 10000 == 0 && el->elem_id > 0, "Generated ", el->elem_id, " ids");

        if (tripplets.size() > BATCH_SIZE * 9 / 10) {
            matrix.insertFromTriplets(tripplets.begin(), tripplets.end());
            logging::info(true, "Squashing buffer");

            tripplets.clear();
        }
    }
    matrix.insertFromTriplets(tripplets.begin(), tripplets.end());
    tripplets.clear();
    return matrix;
}

DynamicVector Model::build_implicit_load_vector(const SparseMatrix &stiffness, const DynamicVector &support) {
    // Size of the support vector
    DynamicVector new_load = DynamicVector::Zero(support.size());

    for (int i = 0; i < support.size(); ++i) {
        // Check if the entry is not NaN
        if (!std::isnan(support[i])) {
            // Loop over the non-zero values of the sparse row
            for (SparseMatrix::InnerIterator it(stiffness, i); it; ++it) {
                new_load[it.row()] += it.value() * support[i];
            }
        }
    }

    return new_load;
}

IndexVector Model::build_mapping_vector(IndexMatrix &unconstrained, IndexMatrix &constrained) {
    // Count non-negative numbers in matrix_a
    int count = (unconstrained.array() >= 0).count();

    IndexVector result(count);
    result.setConstant(-1);  // Initialize with -1 values

    for (int m = 0; m < unconstrained.rows(); ++m) {
        for (int n = 0; n < unconstrained.cols(); ++n) {
            if (unconstrained(m, n) >= 0 && constrained(m, n) >= 0) {
                result(unconstrained(m, n)) = constrained(m, n);
            }
        }
    }

    return result;
}


SparseMatrix Model::build_reduced_stiffness(IndexVector &mapping, SparseMatrix &stiffness) {
    // Find out the new size of the matrix
    int count = (mapping.array() >= 0).count();

    // Use a list of triplets for efficient insertion
    std::vector <Eigen::Triplet<Precision>> tripletList;
    tripletList.reserve(stiffness.nonZeros());

    for (int k = 0; k < stiffness.outerSize(); ++k) {
        for (SparseMatrix::InnerIterator it(stiffness, k); it; ++it) {
            int old_row = it.row();
            int old_col = it.col();

            if (mapping(old_row) != -1 && mapping(old_col) != -1) {

                tripletList.emplace_back(mapping(old_row), mapping(old_col), it.value());
            }
        }
    }

    SparseMatrix reduced_stiffness(count, count);
    reduced_stiffness.setFromTriplets(tripletList.begin(), tripletList.end());

    return reduced_stiffness;
}

DynamicVector Model::build_reduced_load(IndexVector &mapping, DynamicVector &constrained) {
    int count = (mapping.array() >= 0).count();
    DynamicVector reduced_load(count);

    for (int i = 0; i < constrained.size(); ++i) {
        if (mapping(i) != -1) {
            reduced_load(mapping(i)) = constrained(i);
        }
    }

    return reduced_load;
}

NodeData Model::build_global_displacement (IndexMatrix& constrained, DynamicVector& result, IndexMatrix& unconstrained, DynamicVector& supports){
    NodeData res = NodeData(constrained.rows(), constrained.cols());
    res.setZero();

    for (int m = 0; m < constrained.rows(); m++) {
        for (int n = 0; n < constrained.cols(); n++) {
            if (constrained(m,n) >= 0) {
                res(m,n) = result(constrained(m,n));
            }
            if (unconstrained(m,n) >= 0 && !std::isnan(supports(unconstrained(m,n)))){
                res(m,n) = supports(unconstrained(m,n));
            }
        }
    }

    return res;
}


}
}