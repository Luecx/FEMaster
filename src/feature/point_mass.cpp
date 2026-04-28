/**
 * @file point_mass.cpp
 * @brief Implements diagonal point-mass feature assembly.
 *
 * @see src/feature/point_mass.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#include "point_mass.h"

namespace fem {
namespace feature {
namespace {
/**
 * @brief Adds one diagonal matrix entry when the target DOF is active.
 *
 * The global DOF table uses negative ids for inactive DOFs. Zero-valued feature
 * data is also skipped to keep the assembled triplet lists compact.
 *
 * @param global_dof Global DOF id from the system index table.
 * @param value Diagonal value to assemble.
 * @param triplets Triplet list receiving the entry.
 */
void add_diagonal_triplet(int global_dof, Precision value, TripletList& triplets) {
    if (global_dof >= 0 && value != Precision(0)) {
        triplets.emplace_back(global_dof, global_dof, value);
    }
}
} // namespace

/**
 * @copydoc PointMass::assemble_stiffness
 */
void PointMass::assemble_stiffness(const SystemDofIds& indices, TripletList& out) const {
    if (!region_) {
        return;
    }

    for (const ID node_id : *region_) {
        // Translational springs act on ux, uy and uz.
        add_diagonal_triplet(indices(node_id, 0), spring_constants_(0), out);
        add_diagonal_triplet(indices(node_id, 1), spring_constants_(1), out);
        add_diagonal_triplet(indices(node_id, 2), spring_constants_(2), out);

        // Rotational springs act on rx, ry and rz.
        add_diagonal_triplet(indices(node_id, 3), rotary_spring_constants_(0), out);
        add_diagonal_triplet(indices(node_id, 4), rotary_spring_constants_(1), out);
        add_diagonal_triplet(indices(node_id, 5), rotary_spring_constants_(2), out);
    }
}

/**
 * @copydoc PointMass::assemble_mass
 */
void PointMass::assemble_mass(const SystemDofIds& indices, TripletList& out) const {
    if (!region_) {
        return;
    }

    for (const ID node_id : *region_) {
        // Translational point mass uses the same scalar for ux, uy and uz.
        add_diagonal_triplet(indices(node_id, 0), mass_, out);
        add_diagonal_triplet(indices(node_id, 1), mass_, out);
        add_diagonal_triplet(indices(node_id, 2), mass_, out);

        // Rotary inertia is stored component-wise for rx, ry and rz.
        add_diagonal_triplet(indices(node_id, 3), rotary_inertia_(0), out);
        add_diagonal_triplet(indices(node_id, 4), rotary_inertia_(1), out);
        add_diagonal_triplet(indices(node_id, 5), rotary_inertia_(2), out);
    }
}
} // namespace feature
} // namespace fem
