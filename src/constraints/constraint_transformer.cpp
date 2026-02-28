/**
 * @file constraint_transformer.cpp
 * @brief Implements the `ConstraintTransformer` facade utilities.
 *
 * The transformer orchestrates constraint assembly and offers two backends:
 * null-space projection and Lagrange-multiplier (KKT) assembly.
 *
 * @see src/constraints/constraint_transformer.h
 * @see src/constraints/builder/builder.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "constraint_transformer.h"

#include "../core/logging.h"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace fem {
namespace constraint {

namespace {

DynamicVector compute_row_regularization_scale(const SparseMatrix& C) {
    DynamicVector row_scale = DynamicVector::Zero(C.rows());
    for (int col = 0; col < C.outerSize(); ++col) {
        for (SparseMatrix::InnerIterator it(C, col); it; ++it) {
            row_scale[it.row()] += it.value() * it.value();
        }
    }

    for (int i = 0; i < row_scale.size(); ++i) {
        if (!std::isfinite(row_scale[i]) || row_scale[i] <= Precision(0)) {
            row_scale[i] = Precision(1);
        }
    }
    return row_scale;
}

} // namespace

ConstraintTransformer::ConstraintTransformer(const Equations& eqs,
                                             const SystemDofIds& dofs,
                                             Index n_dofs)
    : ConstraintTransformer(eqs, dofs, n_dofs, BuildOptions{}) {}

ConstraintTransformer::ConstraintTransformer(const Equations& eqs,
                                             const SystemDofIds& dofs,
                                             Index n_dofs,
                                             const BuildOptions& opt)
    : method_(opt.method),
      lagrange_regularization_rel_(std::max<Precision>(0, opt.lagrange_regularization_rel)) {
    set_.equations = eqs;
    set_.opt = opt.set;
    set_.assemble(dofs, n_dofs);

    if (method_ == Method::NullSpace) {
        std::tie(map_, report_) = ConstraintBuilder::build(set_, opt.builder);
        return;
    }

    cached_lagrange_lambda_valid_ = false;
    cached_lagrange_lambda_.resize(0);
    lagrange_regularization_row_scale_ = compute_row_regularization_scale(set_.C);
    initialize_identity_map(set_.n);
    report_ = build_lagrange_report();
}

void ConstraintTransformer::initialize_identity_map(Index n) {
    map_.n_ = n;
    map_.r_ = 0;
    map_.nm_ = n;
    map_.masters_.resize(static_cast<std::size_t>(n));
    for (Index i = 0; i < n; ++i) {
        map_.masters_[static_cast<std::size_t>(i)] = i;
    }
    map_.slaves_.clear();
    map_.u_p_ = DynamicVector::Zero(n);
    map_.X_.resize(0, n);

    TripletList trips;
    trips.reserve(static_cast<std::size_t>(n));
    for (Index i = 0; i < n; ++i) {
        trips.emplace_back(i, i, Precision(1));
    }
    map_.T_.resize(n, n);
    map_.T_.setFromTriplets(trips.begin(), trips.end());
    map_.T_.makeCompressed();
}

ConstraintBuilder::Report ConstraintTransformer::build_lagrange_report() const {
    ConstraintBuilder::Report rep;
    rep.m = set_.m;
    rep.n = set_.n;
    rep.homogeneous = (set_.d.size() == 0) || (set_.d.lpNorm<Eigen::Infinity>() == Precision(0));
    rep.feasible = true;
    rep.d_norm = (set_.d.size() == 0) ? Precision(0) : set_.d.norm();
    rep.residual_norm = 0;

    rep.master_idx.resize(static_cast<std::size_t>(set_.n));
    for (Index i = 0; i < set_.n; ++i) {
        rep.master_idx[static_cast<std::size_t>(i)] = i;
    }

    if (set_.C.rows() == 0 || set_.C.cols() == 0) {
        rep.rank = 0;
        rep.n_redundant_rows = 0;
        rep.R11_max_diag = 0;
        return rep;
    }

    // In LAGRANGE mode, skip expensive rank factorization. It is only diagnostic.
    rep.rank = std::min<Index>(rep.m, rep.n);
    rep.n_redundant_rows = std::max<Index>(0, rep.m - rep.rank);
    rep.R11_max_diag = 0;
    return rep;
}

Precision ConstraintTransformer::lagrange_regularization_abs(const SparseMatrix& K) const {
    if (method_ != Method::Lagrange || set_.m == 0 || lagrange_regularization_rel_ <= Precision(0)) {
        return Precision(0);
    }

    Precision sum_abs_diag = 0;
    Index n_diag = 0;
    for (int col = 0; col < K.outerSize(); ++col) {
        for (SparseMatrix::InnerIterator it(K, col); it; ++it) {
            if (it.row() == col) {
                sum_abs_diag += std::abs(it.value());
                ++n_diag;
            }
        }
    }

    Precision k_mean = (n_diag > 0) ? (sum_abs_diag / static_cast<Precision>(n_diag)) : Precision(1);
    if (!std::isfinite(k_mean) || k_mean <= Precision(0)) {
        k_mean = Precision(1);
    }

    return lagrange_regularization_rel_ * k_mean;
}

SparseMatrix ConstraintTransformer::assemble_A(const SparseMatrix& K) const {
    if (method_ == Method::NullSpace) {
        return map_.assemble_A(K);
    }

    logging::error(K.rows() == set_.n && K.cols() == set_.n,
                   "[ConstraintTransformer] K size mismatch for LAGRANGE assemble_A");

    const Index n = set_.n;
    const Index m = set_.m;
    const Index N = n + m;

    TripletList trips;
    trips.reserve(static_cast<std::size_t>(K.nonZeros() + 2 * set_.C.nonZeros() + m));

    for (int col = 0; col < K.outerSize(); ++col) {
        for (SparseMatrix::InnerIterator it(K, col); it; ++it) {
            trips.emplace_back(it.row(), col, it.value());
        }
    }

    for (int col = 0; col < set_.C.outerSize(); ++col) {
        for (SparseMatrix::InnerIterator it(set_.C, col); it; ++it) {
            const Index i = it.row();
            const Index j = col;
            const Precision v = it.value();
            trips.emplace_back(j, n + i, v); // C^T
            trips.emplace_back(n + i, j, v); // C
        }
    }

    const Precision eps_abs = lagrange_regularization_abs(K);
    if (eps_abs > Precision(0)) {
        const bool has_scaled_diag = lagrange_regularization_row_scale_.size() == m;
        for (Index i = 0; i < m; ++i) {
            const Precision scale = has_scaled_diag
                ? lagrange_regularization_row_scale_[i]
                : Precision(1);
            trips.emplace_back(n + i, n + i, -eps_abs * scale);
        }
    }

    SparseMatrix A(N, N);
    A.setFromTriplets(trips.begin(), trips.end());
    A.makeCompressed();
    return A;
}

SparseMatrix ConstraintTransformer::assemble_B(const SparseMatrix& Kg) const {
    if (method_ == Method::NullSpace) {
        return map_.assemble_B(Kg);
    }

    logging::error(false,
                   "[ConstraintTransformer] assemble_B is only available in NULLSPACE mode");
    return SparseMatrix{};
}

DynamicVector ConstraintTransformer::assemble_b(const SparseMatrix& K, const DynamicVector& f) const {
    if (method_ == Method::NullSpace) {
        return map_.assemble_b(K, f);
    }

    logging::error(f.size() == set_.n,
                   "[ConstraintTransformer] f size mismatch for LAGRANGE assemble_b");

    const Index n = set_.n;
    const Index m = set_.m;
    DynamicVector b = DynamicVector::Zero(n + m);
    if (n > 0) {
        b.head(n) = f;
    }
    if (m > 0) {
        logging::error(set_.d.size() == m,
                       "[ConstraintTransformer] d size mismatch while assembling LAGRANGE RHS");
        b.tail(m) = set_.d;
    }
    return b;
}

ConstraintMap::OpA ConstraintTransformer::opA(const SparseMatrix& K) const {
    logging::error(method_ == Method::NullSpace,
                   "[ConstraintTransformer] opA is only available in NULLSPACE mode");
    return ConstraintMap::OpA(map_, K);
}

ConstraintMap::OpB ConstraintTransformer::opB(const SparseMatrix& Kg) const {
    logging::error(method_ == Method::NullSpace,
                   "[ConstraintTransformer] opB is only available in NULLSPACE mode");
    return ConstraintMap::OpB(map_, Kg);
}

void ConstraintTransformer::apply_T(const DynamicVector& q, DynamicVector& u) const {
    if (method_ == Method::NullSpace) {
        map_.apply_T(q, u);
        return;
    }

    u = extract_u_from_solution(q);
}

void ConstraintTransformer::apply_Tt(const DynamicVector& y, DynamicVector& z) const {
    if (method_ == Method::NullSpace) {
        map_.apply_Tt(y, z);
        return;
    }

    z = extract_u_from_solution(y);
}

DynamicVector ConstraintTransformer::extract_u_from_solution(const DynamicVector& x) const {
    const Index n = set_.n;
    const Index m = set_.m;

    if (x.size() == n) {
        return x;
    }

    if (method_ == Method::Lagrange && x.size() == n + m) {
        return x.head(n);
    }

    logging::error(false,
                   "[ConstraintTransformer] Unexpected solution size while extracting displacement");

    DynamicVector u = DynamicVector::Zero(n);
    const Index n_copy = std::min<Index>(n, x.size());
    if (n_copy > 0) {
        u.head(n_copy) = x.head(n_copy);
    }
    return u;
}

DynamicVector ConstraintTransformer::extract_lambda_from_solution(const DynamicVector& x) const {
    const Index n = set_.n;
    const Index m = set_.m;

    if (m == 0) {
        return DynamicVector::Zero(0);
    }

    logging::error(method_ == Method::Lagrange,
                   "[ConstraintTransformer] lambda extraction requested in NULLSPACE mode");

    if (x.size() == n + m) {
        return x.tail(m);
    }

    logging::error(false,
                   "[ConstraintTransformer] Unexpected solution size while extracting multipliers");
    return DynamicVector::Zero(m);
}

void ConstraintTransformer::cache_lagrange_lambda_from_solution(const DynamicVector& x) const {
    cached_lagrange_lambda_valid_ = false;
    cached_lagrange_lambda_.resize(0);
    if (method_ != Method::Lagrange) {
        return;
    }

    if (x.size() == set_.n + set_.m) {
        cached_lagrange_lambda_ = x.tail(set_.m);
        cached_lagrange_lambda_valid_ = true;
    }
}

DynamicVector ConstraintTransformer::recover_u(const DynamicVector& q) const {
    if (method_ == Method::NullSpace) {
        return map_.recover_u(q);
    }
    cache_lagrange_lambda_from_solution(q);
    return extract_u_from_solution(q);
}

DynamicVector ConstraintTransformer::recover_v(const DynamicVector& qdot) const {
    if (method_ == Method::NullSpace) {
        return map_.T() * qdot;
    }
    return extract_u_from_solution(qdot);
}

DynamicVector ConstraintTransformer::recover_a(const DynamicVector& qddot) const {
    if (method_ == Method::NullSpace) {
        return map_.T() * qddot;
    }
    return extract_u_from_solution(qddot);
}

DynamicVector ConstraintTransformer::reactions(const SparseMatrix& K,
                                               const DynamicVector& f,
                                               const DynamicVector& q) const {
    DynamicVector u = recover_u(q);
    return K * u - f;
}

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

    // Tiny diagonal regularization for rank-deficient/near-singular C C^T.
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

DynamicVector ConstraintTransformer::accumulate_support_constraint_forces(const DynamicVector& lambda) const {
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
                g[col] += it.value() * lambda[row];
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
    if (method_ == Method::Lagrange && q.size() == set_.n + set_.m) {
        lambda = extract_lambda_from_solution(q);
    } else {
        lambda = lagrange_multipliers(K, f, q);
    }

    const DynamicVector g = accumulate_support_constraint_forces(lambda);

    // Sign convention for writer output: REACTION_FORCES = K u - f on support DOFs.
    // From equilibrium K u - f + C^T λ = 0, this equals -(C^T λ).
    return -g;
}

void ConstraintTransformer::post_check_static(const SparseMatrix& K,
                                              const DynamicVector& f,
                                              const DynamicVector& u,
                                              Precision tol_constraint_rel,
                                              Precision tol_reduced_rel,
                                              Precision tol_full_rel) const {
    const auto sci6 = [](Precision value) {
        std::ostringstream os;
        os << std::scientific << std::setprecision(6) << value;
        return os.str();
    };

    struct Item {
        std::string label;
        Precision abs_val;
        Precision rel_val;
        Precision tol_val;
        bool pass;
        bool show;
    };

    const auto print_items = [&](const std::vector<Item>& items) {
        int label_width = 22;
        for (const auto& item : items) {
            if (item.show) {
                label_width = std::max<int>(label_width, static_cast<int>(item.label.size()));
            }
        }

        const std::string bullet = "  ";
        const std::string space = " ";
        const std::string colon_sep = " : ";
        const std::string abs_tag = "abs = ";
        const std::string rel_tag = "rel = ";
        const std::string tag_pad(6, ' ');

        logging::info(true, "");
        logging::info(true, "Post-checks");
        logging::up();
        for (const auto& item : items) {
            if (!item.show) {
                continue;
            }

            std::ostringstream line1;
            line1 << bullet
                  << (item.pass ? "[PASS]" : "[FAIL]") << space
                  << std::left << std::setw(label_width) << item.label
                  << colon_sep
                  << abs_tag << sci6(item.abs_val);
            logging::info(true, line1.str());

            std::ostringstream line2;
            line2 << bullet
                  << tag_pad << space
                  << std::left << std::setw(label_width) << ""
                  << colon_sep
                  << rel_tag << sci6(item.rel_val)
                  << "  <= tol = " << sci6(item.tol_val);
            logging::info(true, line2.str());
        }
        logging::down();
    };

    const DynamicVector Cu_d = set().C * u - set().d;
    const Precision abs_Cu_d = Cu_d.norm();
    const Precision den_Cu_d = std::max<Precision>(1, set().d.size() ? set().d.norm() : 0);
    const Precision rel_Cu_d = abs_Cu_d / (den_Cu_d > 0 ? den_Cu_d : 1);

    if (method_ == Method::Lagrange) {
        std::vector<Item> items;
        items.reserve(4);

        items.push_back(Item{"constraints",
                             abs_Cu_d,
                             rel_Cu_d,
                             tol_constraint_rel,
                             rel_Cu_d <= tol_constraint_rel,
                             true});

        const DynamicVector resid = K * u - f;
        if (cached_lagrange_lambda_valid_ && cached_lagrange_lambda_.size() == set_.m) {
            const DynamicVector& lambda = cached_lagrange_lambda_;
            const DynamicVector kkt_u = resid + constraint_forces(lambda);

            const Precision abs_kkt_u = kkt_u.norm();
            const Precision den_kkt_u = std::max<Precision>(1, std::max((K * u).norm(), f.norm()));
            const Precision rel_kkt_u = abs_kkt_u / den_kkt_u;

            items.push_back(Item{"lagrange equilibrium",
                                 abs_kkt_u,
                                 rel_kkt_u,
                                 tol_reduced_rel,
                                 rel_kkt_u <= tol_reduced_rel,
                                 true});

            if (set_.m > 0) {
                const Precision eps_abs = lagrange_regularization_abs(K);
                DynamicVector kkt_c = set_.C * u - set_.d;
                if (eps_abs > Precision(0)) {
                    if (lagrange_regularization_row_scale_.size() == set_.m) {
                        kkt_c.noalias() -= eps_abs * lagrange_regularization_row_scale_.cwiseProduct(lambda);
                    } else {
                        kkt_c.noalias() -= eps_abs * lambda;
                    }
                }

                const Precision abs_kkt_c = kkt_c.norm();
                const Precision den_kkt_c = std::max<Precision>(1, std::max((set_.C * u).norm(), set_.d.norm()));
                const Precision rel_kkt_c = abs_kkt_c / den_kkt_c;

                items.push_back(Item{"lagrange constraints",
                                     abs_kkt_c,
                                     rel_kkt_c,
                                     tol_constraint_rel,
                                     rel_kkt_c <= tol_constraint_rel,
                                     true});
            }
        } else {
            logging::warning(true,
                             "[ConstraintTransformer] Skipping LAGRANGE equilibrium post-check: "
                             "no cached multipliers available.");
        }

        print_items(items);
        return;
    }

    const DynamicVector resid = K * u - f;
    const Precision abs_resid = resid.norm();
    const Precision den_full = std::max<Precision>(1, f.size() ? f.norm() : 0);
    const Precision rel_resid = abs_resid / den_full;

    const Index nm = n_master();
    DynamicVector Ttf(nm);
    DynamicVector TtKu(nm);
    DynamicVector red(nm);
    apply_Tt(f, Ttf);
    apply_Tt(K * u, TtKu);
    apply_Tt(resid, red);

    const Precision abs_red = red.norm();
    const Precision den_red = std::max<Precision>(1, std::max(Ttf.norm(), TtKu.norm()));
    const Precision rel_red = abs_red / den_red;

    std::vector<Item> items;
    items.reserve(5);

    items.push_back(Item{"constraints",
                         abs_Cu_d,
                         rel_Cu_d,
                         tol_constraint_rel,
                         rel_Cu_d <= tol_constraint_rel,
                         true});

    const bool evaluate_full = std::isfinite(tol_full_rel);
    items.push_back(Item{"full residual",
                         abs_resid,
                         rel_resid,
                         tol_full_rel,
                         (!evaluate_full || rel_resid <= tol_full_rel),
                         evaluate_full});

    items.push_back(Item{"reduced equilibrium",
                         abs_red,
                         rel_red,
                         tol_reduced_rel,
                         rel_red <= tol_reduced_rel,
                         true});

    if (nm >= 0) {
        DynamicVector q0 = DynamicVector::Zero(nm);
        DynamicVector u_p = recover_u(q0);
        DynamicVector cup = set().C * u_p - set().d;
        const Precision abs = cup.norm();
        const Precision den = std::max<Precision>(1, set().d.size() ? set().d.norm() : 0);
        const Precision rel = abs / den;
        items.push_back(Item{"affine consistency (u_p)",
                             abs,
                             rel,
                             tol_constraint_rel,
                             rel <= tol_constraint_rel,
                             true});
    }

    if (nm > 0) {
        const int samples = 3;
        Precision abs_acc = 0;
        Precision rel_acc = 0;
        const Precision CnormF = std::max<Precision>(1, set().C.norm());

        for (int s = 0; s < samples; ++s) {
            DynamicVector q = DynamicVector::Random(nm);
            DynamicVector Tq;
            apply_T(q, Tq);
            const Precision un = std::max<Precision>(Precision(1e-16), Tq.norm());
            DynamicVector v = set().C * Tq;
            const Precision an = v.norm();
            abs_acc += an;
            rel_acc += an / std::max<Precision>(1, CnormF * un);
        }

        const Precision abs_avg = abs_acc / samples;
        const Precision rel_avg = rel_acc / samples;
        items.push_back(Item{"null-space (C*T=0)",
                             abs_avg,
                             rel_avg,
                             tol_constraint_rel,
                             rel_avg <= tol_constraint_rel,
                             true});
    }

    print_items(items);
}

} // namespace constraint
} // namespace fem
