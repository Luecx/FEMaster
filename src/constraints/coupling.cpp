/******************************************************************************
 * @file coupling.cpp
 * @brief Implements the Coupling class for handling kinematic _couplings in FEM.
 *
 * This file contains the implementation of the `get_equations` and `master_dofs`
 * methods, which generate the coupling equations and determine the necessary DOFs
 * for master nodes.
 *
 * @see coupling.h
 * @author Finn Eggers
 * @date 04.09.2024
 ******************************************************************************/

#include "coupling.h"

#include "../model/model_data.h"

namespace fem {
namespace constraint {

/******************************************************************************
 * @brief Constructor for the Coupling class.
 *
 * Initializes the coupling constraint with the specified master node, slave nodes,
 * DOFs to couple, and the type of coupling.
 ******************************************************************************/
Coupling::Coupling(ID master_node, model::NodeRegion::Ptr slave_nodes, Dofs coupled_dofs, CouplingType type)
    : master_node(master_node)
    , slave_nodes(slave_nodes)
    , coupled_dofs(coupled_dofs)
    , type(type) {
}

Coupling::Coupling(ID master_node, model::SurfaceRegion::Ptr slave_surfaces, Dofs coupled_dofs, CouplingType type)
    : master_node(master_node)
    , slave_surfaces(slave_surfaces)
    , coupled_dofs(coupled_dofs)
    , type(type) {
}

/******************************************************************************
 * @brief Generates the coupling equations for the specified DOFs.
 *
 * This method computes the coupling equations between the master and slave nodes,
 * handling both translations and rotations. The equations are returned as a list
 * of triplets, which are used to populate the sparse global system matrix.
 *
 * For translational DOFs, the coupling is more complex and involves the rotational
 * DOFs of the master node.
 *
 * @param system_nodal_dofs The global DOF IDs for all nodes in the system.
 * @param node_coords The coordinates of all nodes in the system.
 * @param row_offset The row offset for inserting the equations into the global system.
 * @return TripletList A list of triplets representing the coupling equations.
 ******************************************************************************/
Equations Coupling::get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data) {
    Equations equations {};

    // dont create equations if not kinematic coupling
    if (type != CouplingType::KINEMATIC) {
        return equations;
    }

    auto& node_coords = model_data.get(model::NodeDataEntries::POSITION);

    std::vector<ID> ids{};
    if (slave_nodes != nullptr) {
        for (ID node_id : *slave_nodes) {
            ids.push_back(node_id);
        }
    } else {
        for (ID surface_id : *slave_surfaces) {
            for (ID node_id : *model_data.surfaces[surface_id]) {
                ids.push_back(node_id);
            }
        }
    }

    for (ID slave_node : ids) {
        for (Dim i = 0; i < 6; i++) {
            if (coupled_dofs(i)) {
                if (i >= 3) {
                    // Rotational DOFs: couple with a single value
                    auto dof_master = system_nodal_dofs(master_node, i);
                    auto dof_slave  = system_nodal_dofs(slave_node, i);
                    if (dof_slave < 0)
                        continue;
                    if (dof_master < 0)
                        continue;

                    EquationEntry en1 = EquationEntry {slave_node , i, 1.0};
                    EquationEntry en2 = EquationEntry {master_node, i, -1.0};
                    equations.emplace_back(Equation {en1, en2});
                } else {
                    // Translational DOFs: couple with master translations and rotations
                    auto dof_slave = system_nodal_dofs(slave_node, i);
                    if (dof_slave < 0)
                        continue;

                    Dim u_dof = (Dim) i;
                    Dim r1_dof = (Dim) ((i + 1) % 3 + 3);
                    Dim r2_dof = (Dim) ((i + 2) % 3 + 3);

                    Precision dz            = node_coords(slave_node, 2) - node_coords(master_node, 2);
                    Precision dy            = node_coords(slave_node, 1) - node_coords(master_node, 1);
                    Precision dx            = node_coords(slave_node, 0) - node_coords(master_node, 0);

                    Precision dr1, dr2;
                    if (i == 0) {
                        dr1 = dz;
                        dr2 = -dy;
                    } else if (i == 1) {
                        dr1 = dx;
                        dr2 = -dz;
                    } else {
                        dr1 = dy;
                        dr2 = -dx;
                    }

                    EquationEntry entry_u     = EquationEntry {master_node, u_dof, 1.0};
                    EquationEntry entry_r1    = EquationEntry {master_node, r1_dof, dr1};
                    EquationEntry entry_r2    = EquationEntry {master_node, r2_dof, dr2};
                    EquationEntry entry_slave = EquationEntry {slave_node , u_dof, -1.0};

                    equations.emplace_back(Equation {entry_r1, entry_r2, entry_slave, entry_u});
                }
            }
        }
    }

    return equations;
}


void Coupling::apply_loads(model::ModelData& model_data, NodeData& load_matrix) {
    // Only act for distributing/structural coupling
    if (type != CouplingType::STRUCTURAL) {
        return;
    }

    // Gather slave node IDs (deduplicated)
    std::vector<ID> ids;
    ids.reserve(64);
    if (slave_nodes != nullptr) {
        for (ID node_id : *slave_nodes) ids.push_back(node_id);
    } else if (slave_surfaces != nullptr) {
        for (ID surface_id : *slave_surfaces) {
            for (ID node_id : *model_data.surfaces[surface_id]) {
                ids.push_back(node_id);
            }
        }
    }
    // Deduplicate to avoid over-weighting shared nodes
    {
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    }

    const std::size_t n = ids.size();
    if (n == 0) return;

    // Read master generalized load: [Fx, Fy, Fz, Mx, My, Mz]
    Eigen::Matrix<Precision, 6, 1> b;
    b << load_matrix(master_node, 0), load_matrix(master_node, 1), load_matrix(master_node, 2),
         load_matrix(master_node, 3), load_matrix(master_node, 4), load_matrix(master_node, 5);

    // Nothing to distribute?
    if (b.cwiseAbs().maxCoeff() == Precision(0)) return;

    // Coordinates
    auto& X = model_data.get(model::NodeDataEntries::POSITION);
    const Precision x0 = X(master_node, 0);
    const Precision y0 = X(master_node, 1);
    const Precision z0 = X(master_node, 2);

    // Build A (6 x 3n): [ I I ... ; [r1]x [r2]x ... ], and (optional) weights
    Eigen::Matrix<Precision, Eigen::Dynamic, Eigen::Dynamic> A(6, 3 * (Eigen::Index)n);
    A.setZero();

    // (Uniform) weights: w_i = 1  -> W^{-1} = I   (easy to extend later)
    // If you add weights, prefer scaling A's 3-column blocks by 1/sqrt(w_i) and
    // then build G = As * As^T for improved conditioning.

    auto skew = [](Precision rx, Precision ry, Precision rz) {
        Eigen::Matrix<Precision, 3, 3> S;
        S <<  0,   -rz,  ry,
              rz,   0,  -rx,
             -ry,  rx,   0;
        return S;
    };

    for (std::size_t k = 0; k < n; ++k) {
        const ID i = ids[k];
        const Precision rx = X(i, 0) - x0;
        const Precision ry = X(i, 1) - y0;
        const Precision rz = X(i, 2) - z0;

        // Top block: ... I3 ...
        A.block<3,3>(0, 3 * (Eigen::Index)k).setIdentity();
        // Bottom block: ... [r_i]_x ...
        A.block<3,3>(3, 3 * (Eigen::Index)k) = skew(rx, ry, rz);
    }

    // Build Gram matrix G = A * W^{-1} * A^T; here W^{-1}=I -> G = A*A^T (6x6)
    Eigen::Matrix<Precision, 6, 6> G = (A * A.transpose()).template selfadjointView<Eigen::Lower>();
    // Small Tikhonov in case geometry is rank-deficient
    const Precision eps = std::max(Precision(1e-12), Precision(1e-12) * (G.trace() + Precision(1)));
    G.diagonal().array() += eps;

    // Solve G * lambda = b
    Eigen::LDLT<Eigen::Matrix<Precision, 6, 6>> ldlt(G);
    if (ldlt.info() != Eigen::Success) {
        // As a fallback, increase damping a bit and retry once
        Eigen::Matrix<Precision, 6, 6> Gd = G;
        Gd.diagonal().array() += Precision(1e-9);
        ldlt.compute(Gd);
    }
    Eigen::Matrix<Precision, 6, 1> lambda = ldlt.solve(b);

    // f = W^{-1} * A^T * lambda; with W^{-1}=I -> f = A^T * lambda  (size 3n)
    Eigen::Matrix<Precision, Eigen::Dynamic, 1> f = A.transpose() * lambda;

    // Scatter to slave nodes: add only translational components (Fx,Fy,Fz)
    for (std::size_t k = 0; k < n; ++k) {
        const ID i = ids[k];
        const Eigen::Index j = 3 * (Eigen::Index)k;
        load_matrix(i, 0) += f(j + 0);
        load_matrix(i, 1) += f(j + 1);
        load_matrix(i, 2) += f(j + 2);
        // No nodal moments added here: moment resultant is produced via r x f
    }

    // Consume the master load (it is now represented on slaves)
    for (int c = 0; c < 6; ++c) {
        load_matrix(master_node, c) = Precision(0);
    }
}


/******************************************************************************
 * @brief Computes the necessary DOFs for the master node based on the coupling.
 *
 * This method checks which degrees of freedom (DOFs) need to be active for the
 * master node based on the coupling configuration and the active DOFs in the slave nodes.
 *
 * @param system_dof_mask A mask of active DOFs for all nodes in the system.
 * @return Dofs A Dofs object indicating which DOFs are required for the master node.
 ******************************************************************************/
Dofs Coupling::master_dofs(SystemDofs& system_dof_mask, model::ModelData& model_data) {
    Dofs slave_dofs = {false, false, false, false, false, false};

    // no dofs except for kinematic
    if (type != CouplingType::KINEMATIC) {
        return slave_dofs;
    }

    std::vector<ID> ids{};
    if (slave_nodes != nullptr) {
        for (ID node_id : *slave_nodes) {
            ids.push_back(node_id);
        }
    } else {
        for (ID surface_id : *slave_surfaces) {
            for (ID node_id : *model_data.surfaces[surface_id]) {
                ids.push_back(node_id);
            }
        }
    }


    // Determine which DOFs are active for the slave nodes
    for (ID slave_node : ids) {
        for (int i = 0; i < 6; i++) {
            slave_dofs(i) |= system_dof_mask(slave_node, i);
        }
    }

    // Determine necessary master DOFs based on coupling and slave DOFs
    bool x_translation_needed = coupled_dofs(0) && slave_dofs(0);
    bool y_translation_needed = coupled_dofs(1) && slave_dofs(1);
    bool z_translation_needed = coupled_dofs(2) && slave_dofs(2);

    bool x_rotation_needed    = (coupled_dofs(3) && slave_dofs(3)) || (slave_dofs(1) || slave_dofs(2));
    bool y_rotation_needed    = (coupled_dofs(4) && slave_dofs(4)) || (slave_dofs(0) || slave_dofs(2));
    bool z_rotation_needed    = (coupled_dofs(5) && slave_dofs(5)) || (slave_dofs(0) || slave_dofs(1));

    return Dofs {x_translation_needed,
                 y_translation_needed,
                 z_translation_needed,
                 x_rotation_needed,
                 y_rotation_needed,
                 z_rotation_needed};
}

}    // namespace constraint
}    // namespace fem
