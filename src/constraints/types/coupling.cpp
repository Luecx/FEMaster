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
#include <tuple>
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
 * Distributes the generalized master load of a structural coupling to its
 * slave nodes and, if requested, linearizes the distribution.
 *
 * The master load b = [F; M] is replaced by translational slave forces that
 * reproduce its resultant force and moment about the master position. With
 * lever arms r_k = x_k - x_m and
 *
 *     A = [ I       I      ...  I      ]
 *         [ [r_1]x  [r_2]x ...  [r_n]x ],
 *
 * the minimum-norm distribution is
 *
 *     f = A^T lambda,   G lambda = b,   G = A A^T + eps I,
 *
 * so every slave receives f_k = lambda_f + lambda_m x r_k. The small Tikhonov
 * term eps = 1e-12 (tr(A A^T) + 1) keeps rank-deficient slave patterns solvable.
 *
 * The distribution depends on the positions supplied by model_data. If a load
 * tangent is requested, it is transformed consistently with the load:
 *
 * 1. Master rows db/du of the incoming tangent are consumed together with b
 *    and redistributed through R = A^T G^-1 onto the slave rows.
 * 2. The derivative of the distribution itself with respect to the slave and
 *    master positions is added. For a slave position x_k,
 *
 *        dG/dx_k lambda = H_k,
 *        H_k = [ [lambda_m]x                                      ]
 *              [ -[f_k]x + [r_k]x [lambda_m]x                     ]
 *              + 4e-12 lambda r_k^T,
 *        df_j/dx_k = A_j^T (-G^-1 H_k) + delta_jk [lambda_m]x,
 *
 *    with A_j^T = [I, -[r_j]x]. Because only lever arms depend on the master
 *    position, df/dx_m = -sum_k df/dx_k.
 *
 * The geometric term couples all slave nodes of the coupling with each other
 * and is therefore dense in the slave nodes.
 *
 * @param model_data Positions and surfaces defining the slave pattern.
 * @param load Generalized nodal load field modified in place.
 * @param load_tangent Optional load-tangent triplets transformed in place.
 * @param system_dof_ids Global equation numbering of the load tangent.
 */
void Coupling::apply_loads(model::ModelData&   model_data,
                           model::Field&       load,
                           TripletList*        load_tangent,
                           const SystemDofIds& system_dof_ids) {
    // Only distributing couplings transfer loads
    if (type != CouplingType::STRUCTURAL) {
        return;
    }

    logging::error(model_data.positions != nullptr,
        "COUPLING: positions field is not initialized");

    // Gather the slave nodes. Shared surface nodes are deduplicated so they
    // receive a single share of the distributed load.
    std::vector<ID> ids;
    if (slave_nodes != nullptr) {
        for (ID node_id : *slave_nodes) ids.push_back(node_id);
    } else if (slave_surfaces != nullptr) {
        for (ID surface_id : *slave_surfaces) {
            for (ID node_id : *model_data.surfaces[surface_id]) {
                ids.push_back(node_id);
            }
        }
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

    const Index n = static_cast<Index>(ids.size());
    if (n == 0) return;

    // Read the generalized master load b = [F; M]
    Vec6 b;
    for (Dim c = 0; c < 6; ++c) {
        b(c) = load(master_node, c);
    }

    // Extract the master rows db/du of the incoming load tangent. They are
    // consumed together with b and redistributed onto the slave rows.
    std::vector<std::tuple<Dim, int, Precision>> master_rows;

    if (load_tangent != nullptr) {
        TripletList remaining;
        remaining.reserve(load_tangent->size());

        for (const Triplet& entry : *load_tangent) {
            Dim component = 6;
            for (Dim c = 0; c < 6; ++c) {
                if (system_dof_ids(master_node, c) >= 0 && entry.row() == system_dof_ids(master_node, c)) {
                    component = c;
                }
            }

            if (component < 6) {
                master_rows.emplace_back(component, entry.col(), entry.value());
            } else {
                remaining.push_back(entry);
            }
        }

        load_tangent->swap(remaining);
    }

    if (b.isZero() && master_rows.empty()) return;

    // Lever arms and equilibrium operator A in the supplied configuration
    const auto& positions  = *model_data.positions;
    const Vec3  master_pos = positions.row_vec3(master_node);

    std::vector<Vec3> lever(static_cast<std::size_t>(n));
    DynamicMatrix     A = DynamicMatrix::Zero(6, 3 * n);

    for (Index k = 0; k < n; ++k) {
        lever[k] = positions.row_vec3(ids[k]) - master_pos;

        A.block<3, 3>(0, 3 * k).setIdentity();
        A.block<3, 3>(3, 3 * k) = math::skew(lever[k]);
    }

    // Regularized Gram matrix G = A A^T + eps I. The regularization scales
    // with the trace so that it remains negligible for well-posed patterns.
    StaticMatrix<6, 6> G = A * A.transpose();

    const Precision eps_trace = Precision(1e-12) * (G.trace() + Precision(1));
    const Precision eps       = std::max(Precision(1e-12), eps_trace);
    G.diagonal().array() += eps;

    Eigen::LDLT<StaticMatrix<6, 6>> ldlt(G);
    if (ldlt.info() != Eigen::Success) {
        // Retry once with additional constant damping
        G.diagonal().array() += Precision(1e-9);
        ldlt.compute(G);
    }

    // Minimum-norm slave forces f = A^T G^-1 b
    const Vec6          lambda = ldlt.solve(b);
    const DynamicVector f      = A.transpose() * lambda;

    for (Index k = 0; k < n; ++k) {
        for (Dim d = 0; d < 3; ++d) {
            load(ids[k], d) += f(3 * k + d);
        }
    }

    // Consume the master load; it is now represented on the slaves
    for (Dim c = 0; c < 6; ++c) {
        load(master_node, c) = Precision(0);
    }

    if (load_tangent == nullptr) return;

    // Redistribute the master rows db/du through R = A^T G^-1
    if (!master_rows.empty()) {
        const StaticMatrix<6, 6> G_inv = ldlt.solve(StaticMatrix<6, 6>::Identity());
        const DynamicMatrix      R     = A.transpose() * G_inv;

        for (const auto& [component, col, value] : master_rows) {
            for (Index k = 0; k < n; ++k) {
                for (Dim d = 0; d < 3; ++d) {
                    const int       row          = system_dof_ids(ids[k], d);
                    const Precision contribution = R(3 * k + d, component) * value;
                    if (row < 0 || contribution == Precision(0)) continue;

                    load_tangent->emplace_back(row, col, contribution);
                }
            }
        }
    }

    // The distribution itself does not change if there is no load to distribute
    if (b.isZero()) return;

    // Assemble H_k = dG/dx_k lambda for every slave position. The trace-scaled
    // regularization contributes d(eps)/dx_k = 4e-12 r_k^T.
    const Vec3 lambda_m = lambda.tail<3>();
    const Mat3 skew_m   = math::skew(lambda_m);
    const bool eps_var  = eps_trace >= Precision(1e-12);

    DynamicMatrix H(6, 3 * n);

    for (Index k = 0; k < n; ++k) {
        const Vec3 f_k = f.segment<3>(3 * k);

        H.block<3, 3>(0, 3 * k) = skew_m;
        H.block<3, 3>(3, 3 * k) = -math::skew(f_k) + math::skew(lever[k]) * skew_m;

        if (eps_var) {
            H.block<6, 3>(0, 3 * k) += Precision(4e-12) * lambda * lever[k].transpose();
        }
    }

    // Slave-position derivative df/dx = A^T (-G^-1 H) + blockdiag([lambda_m]x)
    const DynamicMatrix G_inv_H = ldlt.solve(H);
    DynamicMatrix       df_dx   = -(A.transpose() * G_inv_H);

    for (Index k = 0; k < n; ++k) {
        df_dx.block<3, 3>(3 * k, 3 * k) += skew_m;
    }

    // Scatter the slave-slave blocks and the master columns, which follow
    // from df/dx_m = -sum_k df/dx_k because only lever arms depend on x_m
    for (Index j = 0; j < n; ++j) {
        for (Dim d = 0; d < 3; ++d) {
            const int row = system_dof_ids(ids[j], d);
            if (row < 0) continue;

            Vec3 master_column = Vec3::Zero();

            for (Index k = 0; k < n; ++k) {
                for (Dim e = 0; e < 3; ++e) {
                    const Precision value = df_dx(3 * j + d, 3 * k + e);
                    master_column(e) -= value;

                    const int col = system_dof_ids(ids[k], e);
                    if (col < 0 || value == Precision(0)) continue;

                    load_tangent->emplace_back(row, col, value);
                }
            }

            for (Dim e = 0; e < 3; ++e) {
                const int col = system_dof_ids(master_node, e);
                if (col < 0 || master_column(e) == Precision(0)) continue;

                load_tangent->emplace_back(row, col, master_column(e));
            }
        }
    }
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
