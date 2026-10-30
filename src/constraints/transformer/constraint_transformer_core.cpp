/**
 * @file constraint_transformer_core.cpp
 * @brief Implements the method-independent `ConstraintTransformer` interface.
 *
 * This translation unit owns construction of the selected constraint backend
 * and dispatches common matrix, vector and recovery operations to either the
 * affine reduced map or the Lagrange multiplier helpers.
 *
 * @see src/constraints/transformer/constraint_transformer.h
 * @see src/constraints/transformer/constraint_transformer_lagrange.cpp
 * @see src/constraints/transformer/constraint_transformer_postcheck.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "constraint_transformer.h"

#include "../../core/logging.h"
#include "elimination.h"

#include <algorithm>
#include <tuple>
#include <utility>

namespace fem {
namespace constraint {

ConstraintTransformer::ConstraintTransformer(const Equations&    equations,
                                             const SystemDofIds& dofs,
                                             Index               n_dofs)
    : ConstraintTransformer(equations, dofs, n_dofs, Options{}) {}

ConstraintTransformer::ConstraintTransformer(const Equations&    equations,
                                             const SystemDofIds& dofs,
                                             Index               n_dofs,
                                             const Options&      options)
    : method_(options.method),
      lagrange_regularization_(std::max(Precision(0), options.lagrange_regularization)),
      system_(assemble_constraint_system(
          equations,
          dofs,
          n_dofs,
          options.constraints
      )) {
    // Build exactly one backend representation for the assembled system C u = d.
    // Reduced methods store an affine map u = u_p + T q, while the Lagrange
    // method keeps the full displacement space and initializes multiplier data.
    if (method_ == Method::NullSpace) {
        std::tie(map_, report_) = build_null_space(system_, options.null_space);
    } else if (method_ == Method::Elimination) {
        std::tie(map_, report_) = build_elimination(system_);
    } else {
        initialize_lagrange();
    }
}

void ConstraintTransformer::initialize_dof_locations(const SystemDofIds& dofs) {
    logging::error(has_reduced_map(),
        "[ConstraintTransformer] DOF locations require a reduced map");
    logging::error(dofs.cols() > 0,
        "[ConstraintTransformer] DOF numbering has no components");

    // Invert only the master part of the global numbering. This temporary
    // lookup is released when initialization finishes, so the transformer
    // retains no redundant full-system mapping.
    const Index n_master = map_.n_master();
    std::vector<Index> full_to_master(map_.full_size, n_master);

    for (Index master = 0; master < n_master; ++master) {
        const Index gid = map_.masters[master];
        logging::error(gid < map_.full_size && full_to_master[gid] == n_master,
            "[ConstraintTransformer] invalid or repeated master DOF: ", gid);
        full_to_master[gid] = master;
    }

    std::vector<DofLocation> locations(n_master);
    std::vector<Index>       active_nodes;
    std::vector<char>        found(n_master, false);
    active_nodes.reserve(static_cast<std::size_t>(dofs.rows()));

    // Walk the original node-by-component numbering once. A node belongs to
    // the geometric extent even when its active DOFs are all slaves.
    for (Index node = 0; node < static_cast<Index>(dofs.rows()); ++node) {
        bool active = false;
        for (Eigen::Index component = 0; component < dofs.cols(); ++component) {
            const auto gid = dofs(node, component);
            if (gid < 0) continue;

            active = true;
            logging::error(static_cast<Index>(gid) < map_.full_size,
                "[ConstraintTransformer] active DOF index out of range: ", gid);

            const Index master = full_to_master[gid];
            if (master == n_master) continue;

            logging::error(!found[master],
                "[ConstraintTransformer] master DOF appears more than once: ", gid);
            locations[master] = {node, static_cast<int>(component)};
            found[master]     = true;
        }
        if (active) active_nodes.push_back(node);
    }

    // Every reduced coordinate must correspond to one physical nodal DOF.
    for (Index master = 0; master < n_master; ++master) {
        logging::error(found[master],
            "[ConstraintTransformer] no node/DOF for master ", master);
    }

    // Publish the new data only after all numbering checks pass.
    master_locations_          = std::move(locations);
    active_nodes_              = std::move(active_nodes);
    dof_locations_initialized_ = true;
}

const DofLocation& ConstraintTransformer::master_location(Index master) const {
    logging::error(dof_locations_initialized_ &&
                   master < static_cast<Index>(master_locations_.size()),
        "[ConstraintTransformer] master location is unavailable: ", master);
    return master_locations_[master];
}

const std::vector<Index>& ConstraintTransformer::active_nodes() const {
    logging::error(dof_locations_initialized_,
        "[ConstraintTransformer] active nodes require initialize_dof_locations()");
    return active_nodes_;
}

const char* ConstraintTransformer::method_name() const {
    switch (method_) {
        case Method::NullSpace:   return "NULLSPACE";
        case Method::Lagrange:    return "LAGRANGE";
        case Method::Elimination: return "ELIMINATION";
    }

    return "UNKNOWN";
}

const ConstraintMap& ConstraintTransformer::reduced_map() const {
    logging::error(
        has_reduced_map(),
        "[ConstraintTransformer] reduced map requested for LAGRANGE"
    );

    return map_;
}

Index ConstraintTransformer::unknowns() const {
    return method_ == Method::Lagrange
        ? system_.dofs + system_.equations
        : map_.n_master();
}

Index ConstraintTransformer::independent_dofs() const {
    logging::error(
        has_reduced_map(),
        "[ConstraintTransformer] independent DOFs requested for LAGRANGE"
    );

    return map_.n_master();
}

SparseMatrix ConstraintTransformer::assemble_system_matrix(const SparseMatrix& matrix) const {
    // Lagrange produces the augmented saddle-point system. Reduced methods use
    // the congruence transformation T^T K T supplied by the affine map.
    return method_ == Method::Lagrange
        ? assemble_lagrange_matrix(matrix)
        : map_.reduce_matrix(matrix);
}

SparseMatrix ConstraintTransformer::reduce_secondary_matrix(const SparseMatrix& matrix) const {
    logging::error(
        has_reduced_map(),
        "[ConstraintTransformer] secondary matrix reduction is unavailable for LAGRANGE"
    );

    return map_.reduce_matrix(matrix);
}

DynamicVector ConstraintTransformer::assemble_system_rhs(const SparseMatrix&  matrix,
                                                         const DynamicVector& rhs) const {
    // Reduced right-hand sides account for the affine particular displacement
    // through T^T (f - K u_p). Lagrange appends the prescribed constraint values.
    return method_ == Method::Lagrange
        ? assemble_lagrange_rhs(rhs)
        : map_.reduce_rhs(matrix, rhs);
}

void ConstraintTransformer::project_vector(const DynamicVector& full,
                                           DynamicVector&       reduced) const {
    logging::error(
        has_reduced_map(),
        "[ConstraintTransformer] vector projection is unavailable for LAGRANGE"
    );

    map_.project(full, reduced);
}

DynamicVector ConstraintTransformer::project_vector(const DynamicVector& full) const {
    DynamicVector reduced{};
    project_vector(full, reduced);
    return reduced;
}

DynamicVector ConstraintTransformer::recover_displacement(const DynamicVector& solution) const {
    // Absolute displacements include the particular solution for reduced
    // methods. Lagrange solutions already contain the physical displacement
    // block followed by multipliers.
    return method_ == Method::Lagrange
        ? extract_lagrange_displacement(solution)
        : map_.recover(solution);
}

DynamicVector ConstraintTransformer::recover_increment(const DynamicVector& increment) const {
    // Increments, velocities and accelerations are homogeneous quantities, so
    // reduced methods intentionally apply only T and exclude u_p.
    return method_ == Method::Lagrange
        ? extract_lagrange_displacement(increment)
        : map_.T * increment;
}

DynamicVector ConstraintTransformer::recover_velocity(const DynamicVector& velocity) const {
    return recover_increment(velocity);
}

DynamicVector ConstraintTransformer::recover_acceleration(const DynamicVector& acceleration) const {
    return recover_increment(acceleration);
}

DynamicVector ConstraintTransformer::reactions(const SparseMatrix&  matrix,
                                               const DynamicVector& rhs,
                                               const DynamicVector& solution) const {
    return matrix * recover_displacement(solution) - rhs;
}

} // namespace constraint
} // namespace fem
