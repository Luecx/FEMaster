/**
 * @file constraint_map.cpp
 * @brief Implements affine coordinate transformations and reduced-system assembly for constrained systems.
 *
 * `ConstraintMap` represents the affine transformation
 *
 *     u = u_p + T q,
 *
 * where `u` is the full displacement vector, `u_p` is a particular solution
 * satisfying the inhomogeneous constraints, `T` spans the admissible
 * homogeneous displacement space and `q` contains the independent reduced
 * degrees of freedom.
 *
 * The explicit partition into master and slave DOFs allows vector recovery and
 * projection to be performed without multiplying by the complete sparse
 * transformation matrix. Matrix reduction extracts the principal master block
 * directly for selector-only maps and otherwise constructs the general
 * congruence transformation
 *
 *     K_r = T^T K T.
 *
 * @see constraint_map.h
 * @see null_space.cpp
 * @see elimination.cpp
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "constraint_map.h"

#include "../../core/logging.h"

namespace fem {
namespace constraint {

// Recover the complete affine full-system vector
//
//     u = u_p + T q
//
// from a reduced vector q. The explicit master/slave partition is used instead
// of a general sparse matrix-vector multiplication:
//
// - every master DOF receives its corresponding reduced value directly,
// - every slave DOF receives the linear combination stored in X,
// - the particular solution contributes the inhomogeneous affine offset.
void ConstraintMap::recover(const DynamicVector& reduced, DynamicVector& full) const {
    // Verify that the supplied reduced vector contains one entry for every
    // independent master DOF represented by the map.
    logging::warning(
        reduced.size() == static_cast<Eigen::Index>(n_master()),
        "[ConstraintMap::recover] reduced.size()=", reduced.size(),
        " n_master=", n_master()
    );

    // Start from the particular solution u_p. This already contains all
    // constant displacement contributions required by inhomogeneous
    // constraints. For a homogeneous constraint system, this vector is zero.
    full = particular;

    // Insert the identity part of T. Every reduced coordinate q_i corresponds
    // directly to the full-system DOF stored in masters[i], so its contribution
    // is added without consulting the sparse slave transformation X.
    for (Index master_id = 0; master_id < n_master(); ++master_id) {
        full[masters[static_cast<std::size_t>(master_id)]] += reduced[master_id];
    }

    // Add the dependent slave contributions.
    //
    // X has one row per slave DOF and one column per master DOF:
    //
    //     u_slave = X q.
    //
    // Eigen stores the sparse matrix column-wise, so every outer index selects
    // one reduced master coordinate and the inner iterator visits all slave
    // DOFs depending on that coordinate.
    for (Index master_id = 0; master_id < static_cast<Index>(X.outerSize()); ++master_id) {
        for (SparseMatrix::InnerIterator entry(X, master_id); entry; ++entry) {
            // The row of X is a slave-local index. Convert it through the
            // slave partition to the corresponding full-system DOF.
            const Index slave_id = static_cast<Index>(entry.row());

            // Add X(slave_id, master_id) * q(master_id) to the physical slave
            // DOF. The affine contribution from u_p is already present.
            full[slaves[static_cast<std::size_t>(slave_id)]] +=
                entry.value() * reduced[master_id];
        }
    }
}

// Apply the transpose of the homogeneous transformation,
//
//     q_projected = T^T u,
//
// to a full-system vector. This is a projection in the algebraic sense used for
// reduced force and residual vectors; it is not generally the inverse of
// recover(), because T does not have to be orthonormal.
void ConstraintMap::project(const DynamicVector& full, DynamicVector& reduced) const {
    // Verify that the supplied vector is defined over the complete full-system
    // DOF space represented by the map.
    logging::warning(
        full.size() == static_cast<Eigen::Index>(full_size),
        "[ConstraintMap::project] full.size()=", full.size(),
        " full_size=", full_size
    );

    // Resize the result to the number of independent coordinates and clear all
    // entries before accumulating master and slave contributions.
    reduced.setZero(n_master());

    // Apply the identity block of T^T. Each master row of the full vector
    // contributes directly to the associated reduced coordinate.
    for (Index master_id = 0; master_id < n_master(); ++master_id) {
        reduced[master_id] += full[masters[static_cast<std::size_t>(master_id)]];
    }

    // Apply the transposed slave block X^T. For each reduced master coordinate,
    // accumulate the full-vector entries of all dependent slave DOFs weighted
    // by their transformation coefficients.
    for (Index master_id = 0; master_id < static_cast<Index>(X.outerSize()); ++master_id) {
        Precision contribution = Precision(0);

        // Traverse one column of X. Every entry identifies a slave DOF whose
        // full-space value contributes to the current reduced coordinate.
        for (SparseMatrix::InnerIterator entry(X, master_id); entry; ++entry) {
            const Index slave_id = static_cast<Index>(entry.row());

            // Accumulate X(slave_id, master_id) multiplied by the corresponding
            // full-system vector entry.
            contribution += entry.value()
                          * full[slaves[static_cast<std::size_t>(slave_id)]];
        }

        // Add the complete X^T contribution to the direct master contribution
        // already stored in this reduced coordinate.
        reduced[master_id] += contribution;
    }
}

// Allocate and return the full-system representation of a reduced vector. The
// actual affine recovery is delegated to the output-parameter overload so both
// interfaces use exactly the same implementation.
DynamicVector ConstraintMap::recover(const DynamicVector& reduced) const {
    // Allocate storage for every full-system DOF represented by the map.
    DynamicVector full(full_size);

    // Evaluate u = u_p + T q into the allocated vector.
    recover(reduced, full);

    return full;
}

// Allocate and return the reduced projection T^T u of a full-system vector. The
// output-parameter overload performs the actual projection.
DynamicVector ConstraintMap::project(const DynamicVector& full) const {
    // Allocate one entry for every independent master coordinate.
    DynamicVector reduced(n_master());

    // Apply the transpose transformation to the supplied full vector.
    project(full, reduced);

    return reduced;
}

// Reduce a full-system matrix by the congruence transformation
//
//     K_r = T^T K T.
//
// This preserves the virtual work and quadratic-form representation of the
// original operator in the admissible homogeneous coordinate space.
SparseMatrix ConstraintMap::reduce_matrix(const SparseMatrix& matrix) const {
    // A map without slave-to-master coefficients is a pure selector. For the
    // common fixed-DOF case T only embeds the ordered master coordinates, so
    // T^T K T is exactly the principal submatrix K(masters, masters).
    bool ordered_masters = true;
    for (Index master = 1; master < n_master(); ++master) {
        if (masters[static_cast<std::size_t>(master - 1)] >=
            masters[static_cast<std::size_t>(master)]) {
            ordered_masters = false;
            break;
        }
    }

    if (X.nonZeros() == 0 && ordered_masters) {
        logging::error(
            matrix.rows() == static_cast<Eigen::Index>(full_size) &&
            matrix.cols() == static_cast<Eigen::Index>(full_size),
            "[ConstraintMap::reduce_matrix] matrix size does not match full constraint space"
        );

        // Map full-system indices to reduced coordinates. The full_size sentinel
        // marks constrained/slave rows that must be omitted from the principal block.
        const Index invalid = full_size;
        std::vector<Index> full_to_reduced(static_cast<std::size_t>(full_size), invalid);
        for (Index reduced_dof = 0; reduced_dof < n_master(); ++reduced_dof) {
            const Index full_dof = masters[static_cast<std::size_t>(reduced_dof)];
            logging::error(full_dof < full_size,
                "[ConstraintMap::reduce_matrix] master DOF out of range: ", full_dof);
            full_to_reduced[static_cast<std::size_t>(full_dof)] = reduced_dof;
        }

        // Count retained entries per reduced column so Eigen can allocate the
        // output sparsity exactly before the ordered insertion pass.
        Eigen::VectorXi nonzeros_per_column = Eigen::VectorXi::Zero(
            static_cast<Eigen::Index>(n_master())
        );

        for (Index reduced_col = 0; reduced_col < n_master(); ++reduced_col) {
            const Index full_col = masters[static_cast<std::size_t>(reduced_col)];
            int count = 0;
            for (SparseMatrix::InnerIterator entry(matrix, static_cast<int>(full_col)); entry; ++entry) {
                const Index full_row = static_cast<Index>(entry.row());
                if (full_to_reduced[static_cast<std::size_t>(full_row)] != invalid) {
                    ++count;
                }
            }
            nonzeros_per_column[static_cast<Eigen::Index>(reduced_col)] = count;
        }

        SparseMatrix reduced(
            static_cast<Eigen::Index>(n_master()),
            static_cast<Eigen::Index>(n_master())
        );
        reduced.reserve(nonzeros_per_column);

        // The master list is strictly ordered and every source column stores its
        // row indices in ascending order. The mapped reduced rows are therefore
        // also inserted monotonically within each output column.
        for (Index reduced_col = 0; reduced_col < n_master(); ++reduced_col) {
            const Index full_col = masters[static_cast<std::size_t>(reduced_col)];
            for (SparseMatrix::InnerIterator entry(matrix, static_cast<int>(full_col)); entry; ++entry) {
                const Index reduced_row = full_to_reduced[static_cast<std::size_t>(entry.row())];
                if (reduced_row == invalid) continue;

                reduced.insert(
                    static_cast<Eigen::Index>(reduced_row),
                    static_cast<Eigen::Index>(reduced_col)
                ) = entry.value();
            }
        }

        reduced.makeCompressed();
        return reduced;
    }

    // General coupled constraints require the complete congruence transform.
    // Work with compressed sparse operands because Eigen's sparse products are
    // sensitive to the storage state of the input matrices. The transformation
    // transpose is materialized once so the left product does not rebuild it
    // implicitly during the multiplication.
    SparseMatrix matrix_compressed = matrix;
    SparseMatrix T_compressed      = T;
    matrix_compressed.makeCompressed();
    T_compressed.makeCompressed();

    SparseMatrix T_transpose = T_compressed.transpose();
    T_transpose.makeCompressed();

    // Form the right product first. Keeping the intermediate sparse avoids
    // constructing a dense matrix and allows Eigen to exploit the sparsity of
    // both the full operator and the transformation.
    SparseMatrix transformed = matrix_compressed * T_compressed;
    transformed.makeCompressed();

    // Premultiply by T^T to map the full-system rows into reduced coordinates.
    SparseMatrix reduced     = T_transpose * transformed;

    // Convert the result into Eigen's compressed sparse representation before
    // returning it to subsequent assembly or solver operations.
    reduced.makeCompressed();

    return reduced;
}

// Reduce the right-hand side of an affine constrained system.
//
// Substituting
//
//     u = u_p + T q
//
// into
//
//     K u = f
//
// yields
//
//     T^T K T q = T^T (f - K u_p).
//
// The particular solution must therefore be moved to the right-hand side before
// applying the transpose transformation.
DynamicVector ConstraintMap::reduce_rhs(const SparseMatrix&  matrix,
                                        const DynamicVector& rhs) const {
    // Compute the full-system force generated by the prescribed affine
    // displacement component u_p.
    const DynamicVector particular_force = matrix * particular;

    // Remove that contribution from the original load vector so the reduced
    // unknown q only represents the remaining homogeneous displacement part.
    const DynamicVector effective_rhs     = rhs - particular_force;

    // Apply T^T to obtain the right-hand side associated with the reduced
    // equilibrium equations.
    return project(effective_rhs);
}

} // namespace constraint
} // namespace fem
