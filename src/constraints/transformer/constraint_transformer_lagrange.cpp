/**
 * @file constraint_transformer_lagrange.cpp
 * @brief Implements Lagrange-specific `ConstraintTransformer` utilities.
 *
 * @see src/constraints/transformer/constraint_transformer.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "constraint_transformer.h"

#include "../../core/logging.h"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>

#include <algorithm>
#include <cmath>

namespace fem {
namespace constraint {

DynamicVector ConstraintTransformer::solve_multipliers_from_u(const SparseMatrix& K,
                                                              const DynamicVector& f,
                                                              const DynamicVector& u) const {
    if (set_.C.rows() == 0) {
        return DynamicVector::Zero(0);
    }

    const DynamicVector g = f - (K * u);
    const DynamicVector rhs = set_.C * g;

    SparseMatrix normal = set_.C * set_.C.transpose();
    normal.makeCompressed();

    Precision max_diag = 0;
    for (int col = 0; col < normal.outerSize(); ++col) {
        for (SparseMatrix::InnerIterator it(normal, col); it; ++it) {
            if (it.row() == col) {
                max_diag = std::max(max_diag, std::abs(it.value()));
            }
        }
    }
    if (!std::isfinite(max_diag) || max_diag <= Precision(0)) {
        max_diag = Precision(1);
    }

    const Precision reg_abs = std::max<Precision>(Precision(1e-14), Precision(1e-12) * max_diag);
    TripletList reg_trips;
    reg_trips.reserve(static_cast<std::size_t>(set_.m));
    for (Index i = 0; i < set_.m; ++i) {
        reg_trips.emplace_back(i, i, reg_abs);
    }
    SparseMatrix reg_I(set_.m, set_.m);
    reg_I.setFromTriplets(reg_trips.begin(), reg_trips.end());
    normal += reg_I;
    normal.makeCompressed();

    Eigen::SimplicialLDLT<SparseMatrix> ldlt;
    ldlt.compute(normal);
    if (ldlt.info() == Eigen::Success) {
        DynamicVector lambda = ldlt.solve(rhs);
        if (ldlt.info() == Eigen::Success &&
            lambda.size() == set_.C.rows() &&
            lambda.allFinite()) {
            return lambda;
        }
    }

    const SparseMatrix Ct = set_.C.transpose();
    Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr;
    qr.setPivotThreshold(0.0);
    qr.compute(Ct);
    logging::error(qr.info() == Eigen::Success, "[ConstraintTransformer] QR(C^T) factorization failed");

    DynamicVector lambda = qr.solve(g);
    logging::error(lambda.size() == set_.C.rows(), "[ConstraintTransformer] Unexpected lambda size");
    return lambda;
}

DynamicVector ConstraintTransformer::lagrange_multipliers(const SparseMatrix& K,
                                                          const DynamicVector& f,
                                                          const DynamicVector& q) const {
    if (set_.C.rows() == 0) {
        return DynamicVector::Zero(0);
    }

    if (method_ == Method::Lagrange && q.size() == set_.n + set_.m) {
        return extract_lambda_from_solution(q);
    }

    const DynamicVector u = recover_u(q);
    return solve_multipliers_from_u(K, f, u);
}

DynamicVector ConstraintTransformer::constraint_forces(const DynamicVector& lambda) const {
    logging::error(lambda.size() == set_.C.rows(), "[ConstraintTransformer] lambda size mismatch");
    return set_.C.transpose() * lambda;
}

DynamicVector ConstraintTransformer::accumulate_support_constraint_forces(const DynamicVector& lambda,
                                                                          bool use_scaled_rows) const {
    DynamicVector g = DynamicVector::Zero(set_.n);

    if (set_.C.rows() == 0 || set_.equations.empty()) {
        return g;
    }

    logging::error(lambda.size() == set_.C.rows(),
                   "[ConstraintTransformer] lambda size mismatch for support reaction assembly");
    logging::error(static_cast<Index>(set_.kept_row_ids.size()) == set_.C.rows(),
                   "[ConstraintTransformer] kept_row_ids size mismatch");

    for (int col = 0; col < set_.C.outerSize(); ++col) {
        for (SparseMatrix::InnerIterator it(set_.C, col); it; ++it) {
            const int row = it.row();
            const Index eq_idx = set_.kept_row_ids[static_cast<std::size_t>(row)];
            const auto& eq = set_.equations[static_cast<std::size_t>(eq_idx)];
            if (eq.source == EquationSourceKind::Support) {
                const Precision row_scale =
                    (use_scaled_rows && lagrange_row_l2_scale_.size() == set_.m)
                        ? lagrange_row_l2_scale_[row]
                        : Precision(1);
                g[col] += row_scale * it.value() * lambda[row];
            }
        }
    }

    return g;
}

DynamicVector ConstraintTransformer::support_reactions(const SparseMatrix& K,
                                                       const DynamicVector& f,
                                                       const DynamicVector& q) const {
    if (set_.C.rows() == 0 || set_.equations.empty()) {
        return DynamicVector::Zero(set_.n);
    }

    DynamicVector lambda;
    bool use_scaled_rows = false;
    if (method_ == Method::Lagrange && q.size() == set_.n + set_.m) {
        lambda = extract_lambda_from_solution(q);
        use_scaled_rows = true;
    } else {
        lambda = lagrange_multipliers(K, f, q);
    }

    const DynamicVector g = accumulate_support_constraint_forces(lambda, use_scaled_rows);
    return -g;
}
} // namespace constraint
} // namespace fem
