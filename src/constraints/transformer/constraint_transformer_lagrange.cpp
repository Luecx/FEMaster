/**
 * @file constraint_transformer_lagrange.cpp
 * @brief Implements Lagrange-multiplier system assembly and force recovery.
 *
 * The constraint transformer couples structural equilibrium K*u=f to C*u=d
 * through a scaled saddle-point system. This file handles constraint row
 * scaling, multiplier regularization and recovery of constraint forces.
 * ConstraintSystem owns equation data; matrix tools provide sparse block
 * insertion and diagonal magnitudes independently of constraint semantics.
 *
 * @see ConstraintTransformer
 * @see ConstraintSystem
 * @author Finn Eggers
 * @date 11.10.2026
 */

#include "constraint_transformer.h"

#include "../../core/logging.h"
#include "../../mattools/append_sparse_block.h"
#include "../../mattools/max_abs_diagonal.h"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>

#include <algorithm>
#include <cmath>

namespace fem {
namespace constraint {

void ConstraintTransformer::initialize_lagrange() {
    DynamicVector row_norm_squared = DynamicVector::Zero(system_.equations);
    for (int column = 0; column < system_.C.outerSize(); ++column) {
        for (SparseMatrix::InnerIterator entry(system_.C, column); entry; ++entry) {
            row_norm_squared[entry.row()] += entry.value() * entry.value();
        }
    }

    lagrange_row_scale_ = DynamicVector::Ones(system_.equations);
    for (Index row = 0; row < system_.equations; ++row) {
        const Precision norm_squared = row_norm_squared[row];
        if (std::isfinite(norm_squared) && norm_squared > Precision(0)) {
            lagrange_row_scale_[row] = Precision(1) / std::sqrt(norm_squared);
        }
    }

    report_.equations    = system_.equations;
    report_.dofs         = system_.dofs;
    report_.rhs_norm     = system_.d.size() == 0 ? Precision(0) : system_.d.norm();
    report_.residual_norm = Precision(0);
    report_.homogeneous  = system_.d.size() == 0 ||
                           system_.d.lpNorm<Eigen::Infinity>() == Precision(0);
    report_.feasible     = true;
    report_.rank_known   = false;
}

DynamicVector ConstraintTransformer::scale_lagrange_rows(const DynamicVector& values) const {
    if (values.size() == lagrange_row_scale_.size()) {
        return lagrange_row_scale_.cwiseProduct(values);
    }
    return values;
}

/**
 * Determines the multiplier regularization from the structural diagonal scale.
 *
 * Use the maximum absolute stiffness diagonal with a unit fallback for zero or
 * nonfinite scales. No constraints or a disabled regularization return zero.
 * The matrix and transformer state remain unchanged.
 *
 * @param K Structural stiffness matrix in full-system displacement coordinates.
 * @return Configured regularization factor times the usable diagonal scale.
 */
Precision ConstraintTransformer::lagrange_regularization(const SparseMatrix& K) const {
    if (system_.equations == 0 || lagrange_regularization_ <= Precision(0)) {
        return Precision(0);
    }

    // Measure the diagonal scale; retain the local fallback for unusable values.
    Precision stiffness_scale = mattools::max_abs_diagonal(K);
    if (!std::isfinite(stiffness_scale) || stiffness_scale <= Precision(0)) {
        stiffness_scale = Precision(1);
    }
    return lagrange_regularization_ * stiffness_scale;
}

/**
 * Constructs the scaled saddle-point operator [K,C_s^T;C_s,-r*I].
 *
 * Displacement coordinates precede multiplier coordinates. C_s=D*C uses the
 * stored constraint row scales, and r is the stiffness-scaled multiplier
 * regularization. Validate the structural dimensions, collect all block
 * contributions and return independent compressed sparse storage. Inputs and
 * transformer state remain unchanged; a disabled r omits the diagonal block.
 *
 * @param K Full-system structural stiffness with one row and column per DOF.
 * @return Coupled displacement/multiplier matrix.
 */
SparseMatrix ConstraintTransformer::assemble_lagrange_matrix(const SparseMatrix& K) const {
    logging::error(static_cast<Index>(K.rows()) == system_.dofs && static_cast<Index>(K.cols()) == system_.dofs,
        "[Lagrange] stiffness matrix size mismatch");

    const Index total_size = system_.dofs + system_.equations;
    TripletList entries{};
    entries.reserve(static_cast<std::size_t>(
        K.nonZeros() + 2 * system_.C.nonZeros() + system_.equations
    ));

    // Embed structural stiffness in the displacement block of the saddle-point system.
    mattools::append_sparse_block(entries, K);

    for (int column = 0; column < system_.C.outerSize(); ++column) {
        for (SparseMatrix::InnerIterator entry(system_.C, column); entry; ++entry) {
            const Precision value = lagrange_row_scale_[entry.row()] * entry.value();
            entries.emplace_back(column, system_.dofs + entry.row(), value);
            entries.emplace_back(system_.dofs + entry.row(), column, value);
        }
    }

    const Precision regularization = lagrange_regularization(K);
    if (regularization > Precision(0)) {
        for (Index row = 0; row < system_.equations; ++row) {
            entries.emplace_back(system_.dofs + row, system_.dofs + row, -regularization);
        }
    }

    SparseMatrix matrix(total_size, total_size);
    matrix.setFromTriplets(entries.begin(), entries.end());
    matrix.makeCompressed();
    return matrix;
}

DynamicVector ConstraintTransformer::assemble_lagrange_rhs(const DynamicVector& f) const {
    logging::error(static_cast<Index>(f.size()) == system_.dofs,
        "[Lagrange] load vector size mismatch");

    DynamicVector rhs = DynamicVector::Zero(system_.dofs + system_.equations);
    rhs.head(system_.dofs) = f;
    if (system_.equations > 0) {
        rhs.tail(system_.equations) = scale_lagrange_rows(system_.d);
    }
    return rhs;
}

DynamicVector ConstraintTransformer::extract_lagrange_displacement(const DynamicVector& solution) const {
    logging::error(static_cast<Index>(solution.size()) == system_.dofs + system_.equations,
        "[Lagrange] solution size mismatch");
    return solution.head(system_.dofs);
}

DynamicVector ConstraintTransformer::extract_lagrange_multipliers(const DynamicVector& solution) const {
    logging::error(static_cast<Index>(solution.size()) == system_.dofs + system_.equations,
        "[Lagrange] solution size mismatch");
    return solution.tail(system_.equations);
}

/**
 * Recovers multipliers satisfying C^T*lambda=f-K*u in full-system coordinates.
 *
 * Solve the regularized normal equations (C*C^T+r*I)*lambda=C*(f-K*u)
 * using sparse LDLT. Scale r from the largest absolute normal-matrix diagonal
 * with a unit fallback for an unusable scale. If factorization or finite
 * solution recovery fails, solve C^T*lambda=f-K*u through sparse QR instead.
 * An unconstrained system returns an empty vector. Inputs and transformer
 * state are unchanged.
 *
 * @param K Full-system structural stiffness.
 * @param f Full-system external force vector.
 * @param u Full-system displacement vector.
 * @return Multipliers in the original, unscaled constraint row coordinates.
 */
DynamicVector ConstraintTransformer::solve_multipliers(const SparseMatrix& K,
                                                       const DynamicVector& f,
                                                       const DynamicVector& u) const {
    if (system_.equations == 0) {
        return DynamicVector::Zero(0);
    }

    const DynamicVector equilibrium_rhs = f - K * u;
    const DynamicVector normal_rhs      = system_.C * equilibrium_rhs;

    SparseMatrix normal = system_.C * system_.C.transpose();
    normal.makeCompressed();

    // Measure the diagonal scale; retain the local fallback for unusable values.
    Precision max_diagonal = mattools::max_abs_diagonal(normal);
    if (!std::isfinite(max_diagonal) || max_diagonal <= Precision(0)) {
        max_diagonal = Precision(1);
    }

    const Precision regularization = std::max<Precision>(
        Precision(1e-14),
        Precision(1e-12) * max_diagonal
    );
    for (Index row = 0; row < system_.equations; ++row) {
        normal.coeffRef(row, row) += regularization;
    }
    normal.makeCompressed();

    Eigen::SimplicialLDLT<SparseMatrix> ldlt{};
    ldlt.compute(normal);
    if (ldlt.info() == Eigen::Success) {
        DynamicVector multipliers = ldlt.solve(normal_rhs);
        if (ldlt.info() == Eigen::Success && multipliers.allFinite()) {
            return multipliers;
        }
    }

    const SparseMatrix transpose = system_.C.transpose();
    Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr{};
    qr.setPivotThreshold(Precision(0));
    qr.compute(transpose);
    logging::error(qr.info() == Eigen::Success, "[Lagrange] QR(C^T) factorization failed");
    return qr.solve(equilibrium_rhs);
}

DynamicVector ConstraintTransformer::support_constraint_forces(const DynamicVector& multipliers,
                                                               bool scaled_rows) const {
    logging::error(static_cast<Index>(multipliers.size()) == system_.equations,
        "[Lagrange] multiplier size mismatch");
    logging::error(system_.row_sources.size() == static_cast<std::size_t>(system_.equations),
                   "[Lagrange] constraint source count mismatch");

    DynamicVector forces = DynamicVector::Zero(system_.dofs);
    for (int column = 0; column < system_.C.outerSize(); ++column) {
        for (SparseMatrix::InnerIterator entry(system_.C, column); entry; ++entry) {
            if (system_.row_sources[static_cast<std::size_t>(entry.row())] != EquationSourceKind::Support) {
                continue;
            }

            const Precision row_scale = scaled_rows
                                            ? lagrange_row_scale_[entry.row()]
                                            : Precision(1);
            forces[column] += row_scale * entry.value() * multipliers[entry.row()];
        }
    }
    return forces;
}

DynamicVector ConstraintTransformer::support_reactions(const SparseMatrix& K,
                                                       const DynamicVector& f,
                                                       const DynamicVector& solution) const {
    if (system_.equations == 0) {
        return DynamicVector::Zero(system_.dofs);
    }

    DynamicVector multipliers{};
    bool scaled_rows = false;
    if (method_ == Method::Lagrange) {
        multipliers = extract_lagrange_multipliers(solution);
        scaled_rows = true;
    } else {
        multipliers = solve_multipliers(K, f, recover_displacement(solution));
    }

    return -support_constraint_forces(multipliers, scaled_rows);
}

} // namespace constraint
} // namespace fem
