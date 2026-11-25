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

#include <Eigen/QR>
#include <Eigen/SparseQR>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

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
 * Recover support-connected multipliers from C^T*lambda = f-K*u.
 *
 * Rows of C sharing a nonzero column belong to one connected component.
 * Components without support rows cannot contribute to support reactions,
 * even when equations are redundant. Solve the original regularized objective
 * independently on the remaining components without forming C*C^T.
 */
DynamicVector ConstraintTransformer::solve_multipliers(const SparseMatrix& K,
                                                       const DynamicVector& f,
                                                       const DynamicVector& u) const {
    const Index m = system_.equations;
    DynamicVector multipliers = DynamicVector::Zero(m);
    if (m == 0) {
        return multipliers;
    }
    logging::error(system_.row_sources.size() == static_cast<std::size_t>(m),
        "[Support reactions] constraint source count mismatch");

    const auto is_support = [&](Index row) {
        return system_.row_sources[row] == EquationSourceKind::Support;
    };
    if (std::none_of(system_.row_sources.begin(), system_.row_sources.end(),
                     [](EquationSourceKind source) { return source == EquationSourceKind::Support; })) {
        return multipliers;
    }

    const DynamicVector equilibrium_rhs = f - K * u;
    logging::error(equilibrium_rhs.allFinite(),
        "[Support reactions] nonfinite equilibrium residual");

    // Union-find over constraint rows, with union by size and path compression.
    // Each column induces a clique, but only a spanning star is needed.
    std::vector<Index> parent(m);
    std::vector<Index> size(m, 1);
    std::iota(parent.begin(), parent.end(), Index(0));
    const auto root = [&](Index row) {
        while (parent[row] != row) {
            parent[row] = parent[parent[row]];
            row = parent[row];
        }
        return row;
    };
    const auto unite = [&](Index lhs, Index rhs) {
        lhs = root(lhs);
        rhs = root(rhs);
        if (lhs == rhs) return;
        if (size[lhs] < size[rhs]) std::swap(lhs, rhs);
        parent[rhs] = lhs;
        size[lhs] += size[rhs];
    };

    DynamicVector row_norm_squared = DynamicVector::Zero(m);
    for (Index column = 0; column < system_.C.outerSize(); ++column) {
        Index first = m;
        for (SparseMatrix::InnerIterator entry(system_.C, column); entry; ++entry) {
            if (entry.value() == Precision(0)) continue;
            const Index row = entry.row();
            row_norm_squared[row] += entry.value() * entry.value();
            if (first == m) first = row;
            else unite(first, row);
        }
    }
    logging::error(row_norm_squared.allFinite(),
        "[Support reactions] nonfinite constraint row norm");

    // Preserve the previous global regularization convention, including
    // for disconnected components with very different coefficient scales.
    Precision max_diagonal = row_norm_squared.maxCoeff();
    if (max_diagonal <= Precision(0)) max_diagonal = Precision(1);
    const Precision regularization = std::max<Precision>(
        Precision(1e-14), Precision(1e-12) * max_diagonal);
    const Precision damping = std::sqrt(regularization);

    std::vector<char> has_support(m, false);
    for (Index row = 0; row < m; ++row) {
        parent[row] = root(row);
        if (is_support(row)) has_support[parent[row]] = true;
    }

    struct Component {
        std::vector<Index> rows;
        std::vector<Precision> rhs;
        TripletList entries;
    };
    std::vector<Component> blocks;
    std::vector<Index> block_id(m, m);
    std::vector<Index> local_row(m, m);
    Index component_count = 0;
    Index support_components = 0;
    Index selected_rows = 0;
    Index singleton_count = 0;
    Index largest_component = 0;

    for (Index row = 0; row < m; ++row) {
        if (parent[row] != row) continue;
        ++component_count;
        if (!has_support[row]) continue;
        ++support_components;
        selected_rows += size[row];
        largest_component = std::max(largest_component, size[row]);
        if (size[row] == 1) {
            ++singleton_count;
            continue;
        }
        block_id[row] = static_cast<Index>(blocks.size());
        blocks.emplace_back();
        blocks.back().rows.reserve(size[row]);
    }

    for (Index row = 0; row < m; ++row) {
        const Index id = block_id[parent[row]];
        if (id == m) continue;
        auto& block = blocks[id];
        local_row[row] = static_cast<Index>(block.rows.size());
        block.rows.push_back(row);
    }

    logging::info(true, "[Support reactions] components: ", component_count,
                  ", with supports: ", support_components,
                  ", singleton: ", singleton_count);
    logging::info(true, "[Support reactions] selected rows: ", selected_rows,
                  " / ", m, ", largest selected component: ", largest_component);

    // Extract compact local operators directly from C's column-major storage.
    // Each nonzero column belongs to exactly one component. Singleton support
    // rows use the closed form lambda = C_i*b/(||C_i||^2 + regularization).
    for (Index column = 0; column < system_.C.outerSize(); ++column) {
        Index local_dof = 0;
        bool added_column = false;
        for (SparseMatrix::InnerIterator entry(system_.C, column); entry; ++entry) {
            if (entry.value() == Precision(0)) continue;
            const Index row = entry.row();
            const Index representative = parent[row];
            if (!has_support[representative]) break;
            if (size[representative] == 1) {
                multipliers[row] += entry.value() * equilibrium_rhs[column];
                break;
            }

            auto& block = blocks[block_id[representative]];
            if (!added_column) {
                local_dof = static_cast<Index>(block.rhs.size());
                block.rhs.push_back(equilibrium_rhs[column]);
                added_column = true;
            }
            // Store A = C_local^T, scaled by inverse augmented column norms.
            const Precision scale = Precision(1) /
                std::sqrt(row_norm_squared[row] + regularization);
            block.entries.emplace_back(local_dof, local_row[row], entry.value() * scale);
        }
    }

    for (Index row = 0; row < m; ++row) {
        if (has_support[parent[row]] && size[parent[row]] == 1) {
            multipliers[row] /= row_norm_squared[row] + regularization;
        }
    }

    Index dense_solves = 0;
    Index sparse_solves = 0;
    Index largest_dof_count = 0;
    for (auto& block : blocks) {
        const Index n_local = static_cast<Index>(block.rhs.size());
        const Index m_local = static_cast<Index>(block.rows.size());
        largest_dof_count = std::max(largest_dof_count, n_local);

        // Right scaling both blocks keeps the original Tikhonov objective:
        // min_y ||C_local^T D y-b||^2 + regularization*||D y||^2.
        DynamicVector scale(m_local);
        for (Index i = 0; i < m_local; ++i) {
            scale[i] = Precision(1) /
                std::sqrt(row_norm_squared[block.rows[i]] + regularization);
            block.entries.emplace_back(n_local + i, i, damping * scale[i]);
        }

        DynamicVector rhs = DynamicVector::Zero(n_local + m_local);
        for (Index i = 0; i < n_local; ++i) rhs[i] = block.rhs[i];

        DynamicVector y;
        if (m_local <= 32 && n_local <= 256) {
            // Avoid sparse ordering/factorization overhead for tiny components.
            ++dense_solves;
            DynamicMatrix A = DynamicMatrix::Zero(n_local + m_local, m_local);
            for (const auto& entry : block.entries) {
                A(entry.row(), entry.col()) += entry.value();
            }
            y = A.householderQr().solve(rhs);
        } else {
            ++sparse_solves;
            SparseMatrix A(n_local + m_local, m_local);
            A.setFromTriplets(block.entries.begin(), block.entries.end());
            A.makeCompressed();
            Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr;
            qr.setPivotThreshold(Precision(0));
            qr.compute(A);
            logging::error(qr.info() == Eigen::Success &&
                           static_cast<Index>(qr.rank()) == m_local,
                "[Support reactions] local SparseQR factorization failed");
            y = qr.solve(rhs);
            logging::error(qr.info() == Eigen::Success,
                "[Support reactions] local SparseQR solve failed");
        }
        logging::error(y.allFinite(), "[Support reactions] nonfinite local solution");
        for (Index i = 0; i < m_local; ++i) {
            multipliers[block.rows[i]] = scale[i] * y[i];
        }
        // Release each component's assembly storage after solving it.
        block = Component{};
    }
    logging::error(multipliers.allFinite(), "[Support reactions] nonfinite multipliers");
    logging::info(true, "[Support reactions] dense QR solves: ", dense_solves,
                  ", sparse QR solves: ", sparse_solves,
                  ", largest QR component DOFs: ", largest_dof_count);
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
