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

#include <Eigen/IterativeLinearSolvers>

#include <algorithm>
#include <cmath>
#include <limits>

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
 * Recover multipliers from C^T*lambda=f-K*u using damped LSMR.
 *
 * Right-scale the augmented least-squares operator [C^T; sqrt(r)*I]
 * by its inverse column norms. Scaling both blocks preserves the original
 * regularized objective and multiplier coordinates, including redundant rows.
 * Neither C*C^T nor a sparse factorization is constructed.
 */
DynamicVector ConstraintTransformer::solve_multipliers(const SparseMatrix& K,
                                                       const DynamicVector& f,
                                                       const DynamicVector& u) const {
    if (system_.equations == 0) {
        return DynamicVector::Zero(0);
    }

    const DynamicVector equilibrium_rhs = f - K * u;
    logging::error(equilibrium_rhs.allFinite(),
        "[LSMR] nonfinite equilibrium residual");

    // diag(C*C^T) consists of squared row norms; compute it in O(nnz(C)).
    DynamicVector row_norm_squared = DynamicVector::Zero(system_.equations);
    for (Index column = 0; column < system_.C.outerSize(); ++column) {
        for (SparseMatrix::InnerIterator entry(system_.C, column); entry; ++entry) {
            row_norm_squared[entry.row()] += entry.value() * entry.value();
        }
    }
    Precision max_diagonal = row_norm_squared.maxCoeff();
    if (!std::isfinite(max_diagonal) || max_diagonal <= Precision(0)) {
        max_diagonal = Precision(1);
    }
    const Precision regularization = std::max<Precision>(
        Precision(1e-14), Precision(1e-12) * max_diagonal);

    // A right scaling D must also act on the regularization block:
    // min_y ||C^T D y - b||^2 + r ||D y||^2, lambda = D y.
    // Using scalar LSMR damping after scaling would change the objective.
    const DynamicVector scale =
        (row_norm_squared.array() + regularization).sqrt().inverse().matrix();
    logging::error(scale.allFinite() && (scale.array() > Precision(0)).all(),
        "[LSMR] invalid column scaling");

    const Precision damping = std::sqrt(regularization);
    const Index n = system_.dofs;
    const Index m = system_.equations;
    SparseMatrix augmented(n + m, m);
    {
        TripletList entries;
        entries.reserve(static_cast<std::size_t>(system_.C.nonZeros()) +
                        static_cast<std::size_t>(m));
        for (Index column = 0; column < system_.C.outerSize(); ++column) {
            for (SparseMatrix::InnerIterator entry(system_.C, column); entry; ++entry) {
                entries.emplace_back(column, entry.row(), entry.value() * scale[entry.row()]);
            }
        }
        for (Index row = 0; row < m; ++row) {
            entries.emplace_back(n + row, row, damping * scale[row]);
        }
        augmented.setFromTriplets(entries.begin(), entries.end());
    }
    augmented.makeCompressed();

    DynamicVector rhs = DynamicVector::Zero(n + m);
    rhs.head(n) = equilibrium_rhs;
    Eigen::LSMR<SparseMatrix, Eigen::IdentityPreconditioner> solver;
    const Precision tolerance = std::max<Precision>(
        Precision(1e-8), Precision(10) * std::numeric_limits<Precision>::epsilon());
    solver.setTolerance(tolerance);
    solver.setToleranceA(tolerance);
    solver.setToleranceB(tolerance);
    solver.setMaxIterations(4000);
    solver.compute(augmented);
    logging::error(solver.info() == Eigen::Success, "[LSMR] setup failed");

    const DynamicVector scaled_solution = solver.solve(rhs);
    DynamicVector multipliers = scale.cwiseProduct(scaled_solution);
    logging::info(true, "[LSMR] iterations: ", solver.iterations(),
                  ", estimated error: ", solver.error());
    logging::error(multipliers.allFinite(), "[LSMR] nonfinite multipliers");

    // Check the original, unscaled equations independently of Eigen's estimate.
    const DynamicVector residual = equilibrium_rhs - system_.C.transpose() * multipliers;
    const DynamicVector normal_residual = system_.C * residual - regularization * multipliers;
    const Precision rhs_norm = equilibrium_rhs.norm();
    const Precision normal_rhs_norm = (system_.C * equilibrium_rhs).norm();
    logging::error(residual.allFinite() && normal_residual.allFinite(),
        "[LSMR] nonfinite recovery residual");
    logging::info(true, "[LSMR] row norm range: ", std::sqrt(row_norm_squared.minCoeff()),
                  " to ", std::sqrt(row_norm_squared.maxCoeff()));
    logging::info(true, "[LSMR] regularization: ", regularization);
    logging::info(true, "[LSMR] relative equilibrium residual: ",
                  rhs_norm > Precision(0) ? residual.norm() / rhs_norm : residual.norm());
    logging::info(true, "[LSMR] relative regularized normal residual: ",
                  normal_rhs_norm > Precision(0)
                      ? normal_residual.norm() / normal_rhs_norm : normal_residual.norm());
    logging::error(solver.info() == Eigen::Success,
        "[LSMR] multiplier recovery did not converge");
    return multipliers;
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
