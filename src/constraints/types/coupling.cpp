/**
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
 */

#include "coupling.h"

#include "../../math/vec_util.h"
#include "../../model/model_data.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace fem {
namespace constraint {

/**
 * @brief Constructor for the Coupling class.
 *
 * Initializes the coupling constraint with the specified master node, slave nodes,
 * DOFs to couple, and the type of coupling.
 */
Coupling::Coupling(ID master_node, model::NodeRegion::Ptr slave_nodes, Dofs coupled_dofs, CouplingType type)
    : master_node(master_node)
    , slave_nodes(slave_nodes)
    , coupled_dofs(coupled_dofs)
    , type(type) {
    // no-op: reporting pointers set by Model::add_coupling
}

Coupling::Coupling(ID master_node, model::SurfaceRegion::Ptr slave_surfaces, Dofs coupled_dofs, CouplingType type)
    : master_node(master_node)
    , slave_surfaces(slave_surfaces)
    , coupled_dofs(coupled_dofs)
    , type(type) {
    // no-op: reporting pointers set by Model::add_coupling
}

/**
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
 */
Equations Coupling::get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data) {
    Equations equations {};

    if (type == CouplingType::STRUCTURAL) {
        // The old minimum-norm distribution was f_cloud = A^T (A A^T)^-1 f_ref.
        // The work-conjugate kinematics are q_ref = (A A^T)^-1 A u_cloud.
        // The reference geometry keeps these MPC coefficients constant.
        logging::error(model_data.positions_reference != nullptr,
            "COUPLING: reference positions field not set");
        const auto& reference_positions = *model_data.positions_reference;

        std::vector<ID> ids{};
        if (slave_nodes != nullptr) {
            for (ID node_id : *slave_nodes) ids.push_back(node_id);
        } else if (slave_surfaces != nullptr) {
            for (ID surface_id : *slave_surfaces) {
                for (ID node_id : *model_data.surfaces[surface_id]) ids.push_back(node_id);
            }
        }
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        logging::error(!ids.empty(), "COUPLING: distributing cloud is empty");
        logging::error(std::find(ids.begin(), ids.end(), master_node) == ids.end(),
            "COUPLING: reference node cannot be part of its distributing cloud");

        const Index n = static_cast<Index>(ids.size());
        const Vec3 reference = reference_positions.row_vec3(master_node);
        DynamicMatrix A = DynamicMatrix::Zero(6, 3 * n);
        for (Index k = 0; k < n; ++k) {
            const ID node = ids[static_cast<std::size_t>(k)];
            const Vec3 r = reference_positions.row_vec3(node) - reference;
            const Mat3 skew_r = math::skew(r);
            for (Dim d = 0; d < 3; ++d) {
                if (system_nodal_dofs(node, d) < 0) continue;
                A(d, 3 * k + d) = Precision(1);
                A.block<3, 1>(3, 3 * k + d) = skew_r.col(d);
            }
        }

        // Keep the regularization of the former force distribution.
        StaticMatrix<6, 6> G = A * A.transpose();
        const Precision eps_trace = Precision(1e-12) * (G.trace() + Precision(1));
        const Precision eps = std::max(Precision(1e-12), eps_trace);
        G.diagonal().array() += eps;
        Eigen::LDLT<StaticMatrix<6, 6>> ldlt(G);
        if (ldlt.info() != Eigen::Success) {
            G.diagonal().array() += Precision(1e-9);
            ldlt.compute(G);
        }
        logging::error(ldlt.info() == Eigen::Success,
            "COUPLING: failed to factor distributing geometry matrix");

        const DynamicMatrix B = ldlt.solve(A);
        logging::error(B.allFinite(),
            "COUPLING: distributing constraint coefficients are not finite");
        for (Dim component = 0; component < 6; ++component) {
            if (!coupled_dofs(component)) continue;
            logging::error(system_nodal_dofs(master_node, component) >= 0,
                "COUPLING: selected reference DOF is not active");

            Equation equation;
            equation.entries.emplace_back(EquationEntry{master_node, component, Precision(1)});
            for (Index k = 0; k < n; ++k) {
                const ID node = ids[static_cast<std::size_t>(k)];
                for (Dim d = 0; d < 3; ++d) {
                    if (system_nodal_dofs(node, d) < 0) continue;
                    const Precision coeff = -B(component, 3 * k + d);
                    if (coeff != Precision(0)) {
                        equation.entries.emplace_back(EquationEntry{node, d, coeff});
                    }
                }
            }
            logging::error(equation.entries.size() > 1,
                "COUPLING: selected reference DOF cannot be represented by cloud translations");
            equations.push_back(std::move(equation));
        }
        return equations;
    }

    if (type != CouplingType::KINEMATIC) {
        return equations;
    }

    logging::error(model_data.positions != nullptr, "positions field not set in model data");
    const auto& node_coords = *model_data.positions;

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
    {
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
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

                    const Index slave_idx = static_cast<Index>(slave_node);
                    const Index master_idx = static_cast<Index>(master_node);
                    Precision dz            = node_coords(slave_idx, 2) - node_coords(master_idx, 2);
                    Precision dy            = node_coords(slave_idx, 1) - node_coords(master_idx, 1);
                    Precision dx            = node_coords(slave_idx, 0) - node_coords(master_idx, 0);

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


/**
 * @brief Computes the necessary DOFs for the master node based on the coupling.
 *
 * This method checks which degrees of freedom (DOFs) need to be active for the
 * master node based on the coupling configuration and the active DOFs in the slave nodes.
 *
 * @param system_dof_mask A mask of active DOFs for all nodes in the system.
 * @return Dofs A Dofs object indicating which DOFs are required for the master node.
 */
Dofs Coupling::master_dofs(SystemDofs& system_dof_mask, model::ModelData& model_data) {
    Dofs slave_dofs = {false, false, false, false, false, false};

    if (type == CouplingType::STRUCTURAL) {
        return coupled_dofs;
    }
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

static inline std::string dofs_to_string(const Dofs& m) {
    // Order: Ux, Uy, Uz, Rx, Ry, Rz
    static const char* labels[6] = {"Ux","Uy","Uz","Rx","Ry","Rz"};
    std::ostringstream oss;
    bool first = true;
    for (int i = 0; i < 6; ++i) {
        if (m(i)) {
            if (!first) oss << ", ";
            oss << labels[i];
            first = false;
        }
    }
    if (first) oss << "(none)";
    return oss.str();
}

std::string Coupling::str() const {
    std::ostringstream s;
    // Type
    s << (type == CouplingType::KINEMATIC ? "Kinematic" : "Structural") << ": ";
    // Master
    if (master_region) {
        s << "master=" << master_region->name << " (node " << master_node << ")";
    } else {
        s << "master=node " << master_node;
    }
    // Slave
    s << ", slave=";
    if (!slave_nodes && !slave_surfaces) throw std::runtime_error("Coupling::str(): neither slave_nodes nor slave_surfaces set");
    if (slave_nodes)    s << "NSET "  << slave_nodes->name;
    if (slave_surfaces) s << "SFSET " << slave_surfaces->name;
    // DOFs
    s << ", dofs=[" << dofs_to_string(coupled_dofs) << "]";
    return s.str();
}
}    // namespace constraint
}    // namespace fem
