//
// Created by Finn Eggers on 03.09.23.
//

#include "model.h"
#include "../mattools/numerate_dofs.h"
#include "../mattools/assemble.h"
#include "../mattools/assemble_bc.h"

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
    auto res = mattools::numerate_dofs(mask);
    return res;
}
//
// IndexMatrix Model::build_constrained_index_matrix(DynamicVector &support_vector) {
//
//     // first build a boolean matrix of which dofs are present in the system
//     BooleanMatrix mask{this->max_nodes, 6};
//     mask.fill(false);
//
//     // go through each element and mask the dofs for all the nodes
//     for (auto &e: elements) {
//         if (e != nullptr) {
//             auto dofs = e->dofs();
//             for (ID node_local_id = 0; node_local_id < e->n_nodes(); node_local_id++) {
//                 ID node_id = e->nodes()[node_local_id];
//                 for (ID dof = 0; dof < 6; dof++) {
//                     mask(node_id, dof) |= dofs(dof);
//                 }
//             }
//         }
//     }
//
//     // numerate the entries in the mask
//     IndexMatrix res{this->max_nodes, 6};
//
//     ID c = 0;
//     ID c_unconstrained = 0;
//     for (ID n = 0; n < this->max_nodes; n++) {
//         for (ID d = 0; d < 6; d++) {
//             if (mask(n, d)) {
//                 if (std::isnan(support_vector(c_unconstrained))) {
//                     res(n, d) = c++;
//                 } else {
//                     res(n, d) = -1;
//                 }
//                 c_unconstrained++;
//             } else {
//                 res(n, d) = -1;
//             }
//         }
//     }
//
//     return res;
// }

NodeData Model::build_constraint_matrix(std::vector<std::string> supp_sets) {
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

SparseMatrix Model::build_stiffness_matrix(SystemDofIds &indices, ElementData stiffness_scalar) {
    auto lambda = [&](const ElementPtr &el, Precision* storage) {
        MapMatrix stiff = el->stiffness(node_coords, storage);
        stiff *= (stiffness_scalar.rows() > 0) ? stiffness_scalar(el->elem_id, 0) : 1.0;
        return stiff;
    };
    auto res = mattools::assemble_matrix(elements, indices, lambda);
    return res;
}
//
// DynamicVector Model::build_implicit_load_vector(const SparseMatrix &stiffness, const DynamicVector &support) {
//     // Size of the support vector
//     DynamicVector new_load = DynamicVector::Zero(support.size());
//
//     for (int i = 0; i < support.size(); ++i) {
//         // Check if the entry is not NaN
//         if (!std::isnan(support[i])) {
//             // Loop over the non-zero values of the sparse row
//             for (SparseMatrix::InnerIterator it(stiffness, i); it; ++it) {
//                 new_load[it.row()] += it.value() * support[i];
//             }
//         }
//     }
//
//     return new_load;
// }
//
// IndexVector Model::build_mapping_vector(IndexMatrix &unconstrained, IndexMatrix &constrained) {
//     // Count non-negative numbers in matrix_a
//     int count = (unconstrained.array() >= 0).count();
//
//     IndexVector result(count);
//     result.setConstant(-1);  // Initialize with -1 values
//
//     for (int m = 0; m < unconstrained.rows(); ++m) {
//         for (int n = 0; n < unconstrained.cols(); ++n) {
//             if (unconstrained(m, n) >= 0 && constrained(m, n) >= 0) {
//                 result(unconstrained(m, n)) = constrained(m, n);
//             }
//         }
//     }
//
//     return result;
// }
//
//
// SparseMatrix Model::build_reduced_stiffness(IndexVector &mapping, SparseMatrix &stiffness) {
//     // Find out the new size of the matrix
//     int count = (mapping.array() >= 0).count();
//
//     // Use a list of triplets for efficient insertion
//     std::vector <Eigen::Triplet<Precision>> tripletList;
//     tripletList.reserve(stiffness.nonZeros());
//
//     for (int k = 0; k < stiffness.outerSize(); ++k) {
//         for (SparseMatrix::InnerIterator it(stiffness, k); it; ++it) {
//             int old_row = it.row();
//             int old_col = it.col();
//
//             if (mapping(old_row) != -1 && mapping(old_col) != -1) {
//
//                 tripletList.emplace_back(mapping(old_row), mapping(old_col), it.value());
//             }
//         }
//     }
//
//     SparseMatrix reduced_stiffness(count, count);
//     reduced_stiffness.setFromTriplets(tripletList.begin(), tripletList.end());
//
//     return reduced_stiffness;
// }
//
// DynamicVector Model::build_reduced_load(IndexVector &mapping, DynamicVector &constrained) {
//     int count = (mapping.array() >= 0).count();
//     DynamicVector reduced_load(count);
//
//     for (int i = 0; i < constrained.size(); ++i) {
//         if (mapping(i) != -1) {
//             reduced_load(mapping(i)) = constrained(i);
//         }
//     }
//
//     return reduced_load;
// }
//
// NodeData Model::build_global_displacement (IndexMatrix& constrained, DynamicVector& result, IndexMatrix& unconstrained, DynamicVector& supports){
//     NodeData res = NodeData(constrained.rows(), constrained.cols());
//     res.setZero();
//
//     for (int m = 0; m < constrained.rows(); m++) {
//         for (int n = 0; n < constrained.cols(); n++) {
//             if (constrained(m,n) >= 0) {
//                 res(m,n) = result(constrained(m,n));
//             }
//             if (unconstrained(m,n) >= 0 && !std::isnan(supports(unconstrained(m,n)))){
//                 res(m,n) = supports(unconstrained(m,n));
//             }
//         }
//     }
//
//     return res;
// }


}
}