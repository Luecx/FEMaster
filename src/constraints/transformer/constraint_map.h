/**
 * @file constraint_map.h
 * @brief Declares the affine transformation between reduced and full constraint coordinates.
 *
 * `ConstraintMap` stores the representation
 *
 *     u = u_p + T q,
 *
 * used by reduced constraint formulations. The vector `u` contains all active
 * full-system degrees of freedom, `q` contains only the independent reduced
 * coordinates, `u_p` is a particular solution of the inhomogeneous constraint
 * equations and `T` spans the admissible homogeneous displacement space.
 *
 * The full transformation is partitioned into independent master DOFs and
 * dependent slave DOFs. The sparse matrix `X` stores only the slave-to-master
 * block, which allows vector recovery and projection to be evaluated without
 * multiplying by the complete transformation matrix. The complete matrix `T`
 * remains available for congruence transformations of stiffness, mass and
 * other system matrices.
 *
 * @see constraint_map.cpp
 * @see constraint_system.h
 * @see null_space.cpp
 * @see elimination.cpp
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../../core/types_eig.h"

#include <vector>

namespace fem {
namespace constraint {

/**
 * @brief Stores an affine map from independent reduced coordinates to the full system.
 *
 * The map represents every admissible full-system vector in the form
 *
 *     u = u_p + T q.
 *
 * The full DOFs are divided into independent master DOFs and dependent slave
 * DOFs. Each reduced coordinate corresponds directly to one master DOF. The
 * dependent slave values are reconstructed as linear combinations of the
 * reduced coordinates through the sparse matrix `X`.
 *
 * With the full DOFs ordered conceptually as masters and slaves, the
 * homogeneous transformation has the block form
 *
 *     T = [ I ]
 *         [ X ],
 *
 * although the actual full-system rows retain their original global DOF
 * ordering. The vectors `masters` and `slaves` provide the mapping between this
 * logical partition and the original full-system indices.
 *
 * The affine part `particular` is required for inhomogeneous constraints. It
 * satisfies `C u_p = d`, while the columns of `T` satisfy `C T = 0`.
 */
struct ConstraintMap {
    // Number of active DOFs in the original full-system coordinate space. This
    // is the row count of T and the required size of every recovered full
    // vector.
    Index full_size{};

    // Full-system indices of the independent DOFs. Entry i identifies the
    // physical full-system DOF represented directly by reduced coordinate q_i.
    std::vector<Index> masters{};

    // Full-system indices of the dependent DOFs. Entry i identifies the
    // physical full-system DOF represented by row i of the slave
    // transformation block X.
    std::vector<Index> slaves{};

    // Complete homogeneous transformation from reduced coordinates to the
    // original full-system coordinates. Its dimensions are
    //
    //     full_size x n_master().
    //
    // The matrix contains identity entries at the master rows and the
    // slave-to-master coefficients at the slave rows.
    SparseMatrix T{};

    // Slave-to-master transformation block. Its dimensions are
    //
    //     n_slave() x n_master().
    //
    // Each row defines one dependent slave DOF as a linear combination of the
    // independent reduced coordinates.
    SparseMatrix X{};

    // Particular full-system solution u_p of the inhomogeneous constraint
    // equations. This vector contains the constant affine contribution added
    // during displacement recovery. It is zero for homogeneous constraints.
    DynamicVector particular{};

    // Return the number of independent reduced coordinates. This equals the
    // number of master DOFs and the column count of T and X.
    Index n_master() const { return static_cast<Index>(masters.size()); }

    // Return the number of dependent full-system coordinates. This equals the
    // number of slave DOFs and the row count of X.
    Index n_slave() const { return static_cast<Index>(slaves.size()); }

    // Recover a complete full-system vector from reduced coordinates according
    // to
    //
    //     full = particular + T * reduced.
    //
    // The result is written into the supplied output vector.
    void recover(const DynamicVector& reduced, DynamicVector& full) const;

    // Apply the transpose of the homogeneous transformation,
    //
    //     reduced = T^T * full,
    //
    // to a full-system vector. This operation is used to project forces,
    // residuals and other covariant quantities into reduced coordinates.
    void project(const DynamicVector& full, DynamicVector& reduced) const;

    // Allocate and return the complete affine full-system representation of the
    // supplied reduced coordinate vector.
    DynamicVector recover(const DynamicVector& reduced) const;

    // Allocate and return the transpose projection of a full-system vector into
    // the independent reduced coordinate space.
    DynamicVector project(const DynamicVector& full) const;

    // Apply the congruence transformation
    //
    //     reduced = T^T * matrix * T
    //
    // to a full-system operator. This is used for stiffness, mass, geometric
    // stiffness and other matrices acting on displacement-like coordinates.
    SparseMatrix reduce_matrix(const SparseMatrix& matrix) const;

    // Construct the reduced right-hand side for an affine constrained system.
    // The contribution generated by the particular solution is removed before
    // applying the transpose transformation:
    //
    //     reduced_rhs = T^T * (rhs - matrix * particular).
    DynamicVector reduce_rhs(const SparseMatrix& matrix, const DynamicVector& rhs) const;
};

} // namespace constraint
} // namespace fem
