/**
 * @file coupling.h
 * @brief Defines the Coupling class for handling kinematic _couplings in FEM.
 *
 * The Coupling class provides functionality to couple master and slave nodes
 * based on specified degrees of freedom (DOFs). It supports coupling translations
 * and rotations, and generates the necessary coupling equations for FEM systems.
 *
 * @note Additional coupling types can be added as needed.
 *
 * @see coupling.cpp
 * @author Finn Eggers
 * @date 04.09.2024
 */

#pragma once    // Ensures this file is only included once during compilation

#include "../core/types_eig.h"
#include "../core/types_cls.h"
#include "../data/region.h"
#include "equation.h"
#include "../core/printable.h"
#include <string>

namespace fem {
namespace constraint {

/**
 * @enum CouplingType
 * @brief Defines the available types of _couplings in FEM.
 *
 * Currently, only kinematic coupling is supported, but additional types can be
 * added in the future.
 */
enum class CouplingType {
    KINEMATIC,    ///< Kinematic coupling between master and slave nodes
    STRUCTURAL,   ///< Structural distribution coupling
};

/**
 * @class Coupling
 * @brief Implements a coupling constraint between master and slave nodes.
 *
 * The Coupling class generates the necessary equations to couple the degrees
 * of freedom (DOFs) of slave nodes with a master node, based on a specified
 * coupling type and the system's DOF configurations.
 */
class Coupling : public fem::Printable {

    public:
    ID                        master_node;                 ///< Master node ID
    model::NodeRegion::Ptr    slave_nodes    = nullptr;    ///< List of slave node IDs
    model::SurfaceRegion::Ptr slave_surfaces = nullptr;    ///< List of slave surface IDs
    Dofs                      coupled_dofs;                ///< DOFs that are to be coupled (6 DOFs per node)
    CouplingType              type;                        ///< Type of coupling (e.g., kinematic)

    // Reporting: keep pointer to the original master region for its name
    model::NodeRegion::Ptr    master_region = nullptr;    ///< Name and id of the master set (expects size==1)

    public:
    /**
     * @brief Constructor for the Coupling class.
     *
     * Initializes the coupling constraint with a master node, a set of slave nodes,
     * the degrees of freedom to couple, and the type of coupling.
     *
     * @param master_node The ID of the master node.
     * @param slave_nodes A vector of IDs representing the slave nodes.
     * @param coupled_dofs A Dofs object specifying the DOFs to couple.
     * @param type The type of coupling to apply (e.g., kinematic).
     */
    Coupling(ID master_node, model::NodeRegion::Ptr slave_nodes, Dofs coupled_dofs, CouplingType type);
    Coupling(ID master_node, model::SurfaceRegion::Ptr slave_surfaces, Dofs coupled_dofs, CouplingType type);

    /**
     * @brief Generates the coupling equations for the coupled DOFs.
     *
     * This function generates the necessary equations to couple the slave node
     * DOFs with the master node, based on translations and rotations. The
     * equations are returned as a list of triplets, representing matrix entries
     * in a sparse system.
     *
     * @param system_nodal_dofs The global DOF IDs for all nodes in the system.
     * @param model_data model data which has the coordinates of all nodes in the system.
     * @param row_offset The row offset for inserting the equations into the global system.
     * @return TripletList A list of triplets representing the coupling equations.
     */
    Equations get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data);

    /**
     * @brief Distribute loads from the master to slaves with static equivalence.
     *
     * Reads the generalized load stored at @p load_matrix for @p master_node
     * (translations -> forces, rotations -> moments), removes it (so it is not
     * applied at the master), and scatters **equivalent nodal forces** to the
     * slave nodes so that:
     * \f[
     *   \sum_i f_i = F,\qquad \sum_i (x_i - x_0) \times f_i = M,
     * \f]
     * where \f$x_0\f$ is the master position and \f$x_i\f$ are slave node
     * positions obtained from @p model_data.
     *
     * Typical implementation uses a weighted minimum-norm map
     * \f$ f = W^{-1} A^\top (A W^{-1} A^\top)^{-1} [F;M] \f$ with
     * \f$A=[I\ I\ \dots;\ [r_1]_\times\ [r_2]_\times\ \dots]\f$ and
     * \f$r_i=x_i-x_0\f$, where \f$W\f$ encodes tributary areas/lengths.
     * For surface slaves, a consistent element-wise load integration is a good
     * alternative (constant + linear traction field reproduction).
     *
     * @param model_data  Provides master/slave coordinates and mesh topology.
     * @param load_matrix In/out per-node load storage (RHS): the master's load
     *                    is consumed and equivalent slave nodal loads are added.
     *
     * @pre @p load_matrix contains the master's intended load (Fx,Fy,Fz,Mx,My,Mz).
     * @pre Slave node set derivable from @c slave_nodes or @c slave_surfaces.
     * @post @p load_matrix[master_node] is cleared (or reduced by what was
     *       redistributed); slave nodes receive the distributed forces.
     *
     * @warning If slave geometry is rank-deficient (for example, all slaves on the
     *          master point), some moments cannot be reproduced; implementation
     *          should detect this and either damp/invert in the attainable
     *          subspace or issue a warning.
     * @note This method **does not** add constraints; it only assembles the RHS.
     * @note Use this also when you conceptually "apply a load to a reference
     *       point" but want the structure to carry it through a patch of nodes.
     */
    void apply_loads(model::ModelData& model_data, NodeData& load_matrix);

    /**
     * @brief Computes the necessary DOFs for the master node based on the coupling.
     *
     * This function determines which degrees of freedom (DOFs) are required for the
     * master node, based on the coupled DOFs and the active DOFs in the slave nodes.
     *
     * @param system_dof_mask A mask of active DOFs for all nodes in the system.
     * @return Dofs A Dofs object indicating which DOFs are required for the master node.
     */
    Dofs master_dofs(SystemDofs& system_dof_mask, model::ModelData& model_data);

    /**
     * @brief One-line string representation describing this coupling.
     *
     * Includes coupling type, master (set name or node id), slave set name,
     * whether it is a node or surface set, and the coupled DOFs.
     */
    std::string str() const override;
};

}    // namespace constraint
}    // namespace fem
