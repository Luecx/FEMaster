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
#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

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

    // dont create equations if not kinematic coupling
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
 * @brief Distributes a generalized master load and optionally its tangent.
 *
 * The structural coupling does not introduce constraint equations: it maps
 * the master force/moment to translational slave forces. In nonlinear assembly
 * the derivative must include changes of the incoming master force and the
 * geometry-dependent distribution map itself.
 */
void Coupling::apply_loads(
    model::ModelData& model_data,
    model::Field& load_matrix,
    const SystemDofIds* indices,
    SparseMatrix* K_external
) {
    if (type != CouplingType::STRUCTURAL) return;

    logging::error(K_external == nullptr || indices != nullptr,
        "STRUCTURAL COUPLING: external tangent requires DOF numbering");

    // Collect every slave exactly once, irrespective of how surface sets overlap.
    std::vector<ID> ids;
    if (slave_nodes) {
        for (ID node : *slave_nodes) ids.push_back(node);
    } else if (slave_surfaces) {
        for (ID surface : *slave_surfaces) {
            for (ID node : *model_data.surfaces[surface]) ids.push_back(node);
        }
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    if (ids.empty()) return;

    Vec6 master_load;
    for (Index k = 0; k < 6; ++k) {
        master_load[k] = load_matrix(master_node, k);
    }
    if (!K_external && master_load.isZero()) return;

    logging::error(model_data.positions != nullptr,
        "STRUCTURAL COUPLING: current positions are unavailable");
    const auto& positions = *model_data.positions;
    Vec3 master_position = positions.row_vec3(master_node);

    std::vector<Vec3> slave_positions;
    slave_positions.reserve(ids.size());
    for (ID node : ids) slave_positions.push_back(positions.row_vec3(node));

    const Index n_slaves = static_cast<Index>(ids.size());

    // The least-norm mapping D = A^T (A A^T + eps I)^-1 is identical to
    // the previous force-only implementation. Its geometry arguments are local
    // so differentiation never modifies global model coordinates.
    auto distribution = [&](const Vec3& center, const std::vector<Vec3>& nodes) {
        DynamicMatrix A = DynamicMatrix::Zero(6, 3 * n_slaves);
        for (Index i = 0; i < n_slaves; ++i) {
            const Vec3 r = nodes[static_cast<std::size_t>(i)] - center;
            A.block<3, 3>(0, 3 * i).setIdentity();
            A.block<3, 3>(3, 3 * i) = math::skew(r);
        }

        Mat6 gram = (A * A.transpose()).template selfadjointView<Eigen::Lower>();
        const Precision eps = std::max(
            Precision(1e-12), Precision(1e-12) * (gram.trace() + Precision(1)));
        gram.diagonal().array() += eps;

        Eigen::LDLT<Mat6> ldlt(gram);
        if (ldlt.info() != Eigen::Success) {
            gram.diagonal().array() += Precision(1e-9);
            ldlt.compute(gram);
        }
        logging::error(ldlt.info() == Eigen::Success,
            "STRUCTURAL COUPLING: failed to factor distribution matrix");

        // Solve all six unit load cases to obtain the complete redistribution map.
        return DynamicMatrix(A.transpose() * ldlt.solve(Mat6::Identity()));
    };

    const DynamicMatrix D = distribution(master_position, slave_positions);
    const DynamicVector slave_load = D * master_load;

    // Redistribute the force itself, retaining any other loads already on slaves.
    for (Index i = 0; i < n_slaves; ++i) {
        for (Index d = 0; d < 3; ++d) {
            load_matrix(ids[static_cast<std::size_t>(i)], d) += slave_load[3 * i + d];
        }
    }
    for (Index d = 0; d < 6; ++d) load_matrix(master_node, d) = Precision(0);

    if (!K_external) return;

    // Push forward every incoming master-force derivative:
    //     d(D b)/du = D db/du + (dD/du) b.
    // Other rows, including independently applied slave loads, are unchanged.
    TripletList triplets;
    triplets.reserve(static_cast<std::size_t>(K_external->nonZeros())
                     + static_cast<std::size_t>(n_slaves * 18));
    std::map<Index, Vec6> master_derivatives;

    for (Index outer = 0; outer < K_external->outerSize(); ++outer) {
        for (SparseMatrix::InnerIterator entry(*K_external, outer); entry; ++entry) {
            Index master_component = -1;
            for (Index d = 0; d < 6; ++d) {
                if ((*indices)(master_node, d) >= 0
                    && entry.row() == (*indices)(master_node, d)) {
                    master_component = d;
                    break;
                }
            }
            if (master_component < 0) {
                triplets.emplace_back(entry.row(), entry.col(), entry.value());
                continue;
            }

            auto [it, inserted] =
                master_derivatives.try_emplace(entry.col(), Vec6::Zero());
            it->second[master_component] += entry.value();
        }
    }

    for (const auto& [column, derivative] : master_derivatives) {
        const DynamicVector distributed = D * derivative;
        for (Index i = 0; i < n_slaves; ++i) {
            for (Index d = 0; d < 3; ++d) {
                const Index row = (*indices)(ids[static_cast<std::size_t>(i)], d);
                if (row >= 0 && distributed[3 * i + d] != Precision(0)) {
                    triplets.emplace_back(row, column, distributed[3 * i + d]);
                }
            }
        }
    }

    // Differentiate the distribution with the incoming master load fixed.
    // Only the master/slave translational coordinates affect the local map.
    // A central difference includes the derivative of the regularized Gram
    // solve without constructing or clearing global scratch fields.
    if (!master_load.isZero()) {
        std::vector<ID> participants = ids;
        participants.push_back(master_node);
        std::sort(participants.begin(), participants.end());
        participants.erase(std::unique(participants.begin(), participants.end()),
                           participants.end());

        const Precision step_scale =
            std::cbrt(std::numeric_limits<Precision>::epsilon());

        for (ID node : participants) {
            for (Index d = 0; d < 3; ++d) {
                const Index column = (*indices)(node, d);
                if (column < 0) continue;

                const Precision x = positions(node, d);
                const Precision h = step_scale * std::max(Precision(1), std::abs(x));

                // The master can in principle also occur in the slave region.
                // Perturb every occurrence of the same physical node together.
                if (node == master_node) master_position[d] += h;
                for (Index i = 0; i < n_slaves; ++i) {
                    if (ids[static_cast<std::size_t>(i)] == node)
                        slave_positions[static_cast<std::size_t>(i)][d] += h;
                }
                const DynamicVector plus =
                    distribution(master_position, slave_positions) * master_load;

                if (node == master_node) master_position[d] -= Precision(2) * h;
                for (Index i = 0; i < n_slaves; ++i) {
                    if (ids[static_cast<std::size_t>(i)] == node)
                        slave_positions[static_cast<std::size_t>(i)][d] -= Precision(2) * h;
                }
                const DynamicVector minus =
                    distribution(master_position, slave_positions) * master_load;

                if (node == master_node) master_position[d] += h;
                for (Index i = 0; i < n_slaves; ++i) {
                    if (ids[static_cast<std::size_t>(i)] == node)
                        slave_positions[static_cast<std::size_t>(i)][d] += h;
                }

                const DynamicVector derivative =
                    (plus - minus) / (Precision(2) * h);
                for (Index i = 0; i < n_slaves; ++i) {
                    for (Index k = 0; k < 3; ++k) {
                        const Index row = (*indices)(ids[static_cast<std::size_t>(i)], k);
                        if (row >= 0 && derivative[3 * i + k] != Precision(0))
                            triplets.emplace_back(row, column, derivative[3 * i + k]);
                    }
                }
            }
        }
    }

    K_external->setZero();
    K_external->setFromTriplets(triplets.begin(), triplets.end());
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
