/******************************************************************************
 * @file constraint_transformer.cpp
 * @brief Diagnostics & checklist for static post-checks on ConstraintTransformer.
 ******************************************************************************/

#include "constraint_transformer.h"

#include "../core/logging.h"

#include <cmath>     // std::isfinite
#include <iomanip>   // std::setw, std::setprecision, std::scientific
#include <sstream>   // std::ostringstream
#include <vector>
#include <string>
#include <algorithm>

namespace fem { namespace constraint {

// -----------------------------------------------------------------------------
// Convenience overloads (default options)
// -----------------------------------------------------------------------------
ConstraintTransformer::StaticCheckResult
ConstraintTransformer::check_static(const SparseMatrix& K,
                                    const DynamicVector& f,
                                    const DynamicVector& u) const
{
    return check_static(K, f, u, StaticCheckOptions{});
}

void ConstraintTransformer::print_checklist(const StaticCheckResult& r) const {
    print_checklist(r, StaticCheckOptions{});
}

// -----------------------------------------------------------------------------
// Main implementations
// -----------------------------------------------------------------------------
ConstraintTransformer::StaticCheckResult
ConstraintTransformer::check_static(const SparseMatrix& K,
                                    const DynamicVector& f,
                                    const DynamicVector& u,
                                    StaticCheckOptions opt) const
{
    StaticCheckResult r;

    // 1) Norm of u (scale indicator)
    r.norm_u = u.norm();

    // 2) Constraint satisfaction
    DynamicVector cu_d = this->set().C * u - this->set().d;
    r.norm_Cu_d  = cu_d.norm();
    r.denom_Cu_d = std::max<Precision>(1, this->set().d.size() ? this->set().d.norm() : 0);
    r.rel_Cu_d   = (r.denom_Cu_d > 0) ? (r.norm_Cu_d / r.denom_Cu_d) : r.norm_Cu_d;
    r.pass_constraints = (r.rel_Cu_d <= opt.tol_constraint_rel);

    // 3) Full residual
    DynamicVector resid = K * u - f;
    r.norm_resid = resid.norm();
    r.denom_full = std::max<Precision>(1, f.size() ? f.norm() : 0);
    r.rel_resid  = (r.denom_full > 0) ? (r.norm_resid / r.denom_full) : r.norm_resid;

    // 4) Reduced (projected) residual: z = T^T y
    const Index nm = n_master();
    DynamicVector Ttf(nm), TtKu(nm), red(nm);
    apply_Tt(f,      Ttf);
    apply_Tt(K * u,  TtKu);
    apply_Tt(resid,  red);

    r.norm_reduced  = red.norm();
    r.denom_reduced = std::max<Precision>(1, std::max(Ttf.norm(), TtKu.norm()));
    r.rel_reduced   = (r.denom_reduced > 0) ? (r.norm_reduced / r.denom_reduced) : r.norm_reduced;
    r.pass_reduced_eq = (r.rel_reduced <= opt.tol_reduced_rel);

    // 5) Optional full residual gate
    if (std::isfinite(opt.tol_full_rel))
        r.pass_full_resid = (r.rel_resid <= opt.tol_full_rel);

    return r;
}

void ConstraintTransformer::print_checklist(const StaticCheckResult& r,
                                            StaticCheckOptions opt) const
{
    // number formatting
    auto sci6 = [](Precision x) {
        std::ostringstream os;
        os << std::scientific << std::setprecision(6) << x;
        return os.str();
    };

    struct Item {
        std::string label;
        Precision   abs_val;
        Precision   rel_val;
        bool        show_tol;
        Precision   tol_val;
        bool        pass;
        bool        show; // whether to include in this run
    };

    // Build items in the order we want to print them
    std::vector<Item> items;
    items.push_back(Item{
        "constraints",
        r.norm_Cu_d, r.rel_Cu_d,
        /*show_tol*/true, opt.tol_constraint_rel,
        r.pass_constraints,
        /*show*/true
    });

    items.push_back(Item{
        "full residual",
        r.norm_resid, r.rel_resid,
        /*show_tol*/true, opt.tol_full_rel,
        std::isfinite(opt.tol_full_rel) ? (r.rel_resid <= opt.tol_full_rel) : true,
        /*show*/std::isfinite(opt.tol_full_rel)
    });

    items.push_back(Item{
        "reduced equilibrium",
        r.norm_reduced, r.rel_reduced,
        /*show_tol*/true, opt.tol_reduced_rel,
        r.pass_reduced_eq,
        /*show*/true
    });

    // Extra checks:
    // (A) Affine consistency: C * u_p ≈ d
    {
        const Index nm = n_master();
        if (nm >= 0) {
            DynamicVector q0 = DynamicVector::Zero(nm);
            DynamicVector u_p = recover_u(q0);                 // u_p + T*0 = u_p
            DynamicVector cup = set().C * u_p - set().d;

            const Precision abs  = cup.norm();
            const Precision den  = std::max<Precision>(1, set().d.size() ? set().d.norm() : 0);
            const Precision rel  = (den > 0) ? (abs / den) : abs;
            const bool      pass = (rel <= opt.tol_constraint_rel);

            items.push_back(Item{
                "affine consistency (u_p)",
                abs, rel,
                /*show_tol*/true, opt.tol_constraint_rel,
                pass,
                /*show*/true
            });
        }
    }

    // (B) Null-space property (stochastic): C * (T q) ≈ 0  -> ASCII label: "null-space (C*T=0)"
    {
        const Index nm = n_master();
        if (nm > 0) {
            const int    S      = 3; // samples
            Precision    abs_acc = 0;
            Precision    rel_acc = 0;

            // Frobenius norm proxy from Eigen::SparseMatrix::norm() (compatible with double)
            const Precision CnormF = std::max<Precision>(1, set().C.norm());

            for (int s = 0; s < S; ++s) {
                DynamicVector q = DynamicVector::Random(nm);
                DynamicVector Tq;
                apply_T(q, Tq);
                const Precision un = std::max<Precision>(1e-16, Tq.norm());

                DynamicVector v = set().C * Tq;       // should be ~0
                const Precision an = v.norm();

                abs_acc += an;
                rel_acc += an / std::max<Precision>(1, CnormF * un);
            }

            const Precision abs_avg = abs_acc / S;
            const Precision rel_avg = rel_acc / S;
            const bool      pass    = (rel_avg <= opt.tol_constraint_rel);

            items.push_back(Item{
                "null-space (C*T=0)",
                abs_avg, rel_avg,
                /*show_tol*/true, opt.tol_constraint_rel,
                pass,
                /*show*/true
            });
        }
    }

    // Compute dynamic label width from the items we actually print
    std::size_t max_label_len = 0;
    for (const auto& it : items) {
        if (!it.show) continue;
        max_label_len = std::max<std::size_t>(max_label_len, it.label.size());
    }
    // Enforce a sensible minimum so short labels still look nice
    const int label_w = static_cast<int>(std::max<std::size_t>(22, max_label_len));

    // Constants for printing
    const std::string bullet    = "  ";
    const std::string sp        = " ";
    const std::string colon_sep = " : ";
    const std::string abs_tag   = "abs = ";
    const std::string rel_tag   = "rel = ";
    const std::string tag_pad(6, ' '); // width of "[PASS]" or "[FAIL]"

    // Helper to print a two-line block with a single PASS/FAIL tag on the first line.
    auto print_block = [&](const Item& it)
    {
        // Line 1: "[PASS] Label............ : abs = <num>"
        {
            std::ostringstream os;
            os << bullet
               << (it.pass ? "[PASS]" : "[FAIL]") << sp
               << std::left << std::setw(label_w) << it.label
               << colon_sep
               << abs_tag << sci6(it.abs_val);
            logging::info(true, os.str());
        }

        // Line 2: aligned to the same colon, no PASS/FAIL tag, no extra text
        {
            std::ostringstream os;
            os << bullet
               << tag_pad << sp
               << std::left << std::setw(label_w) << ""    // empty label space
               << colon_sep
               << rel_tag << sci6(it.rel_val);

            if (it.show_tol) {
                os << "  <= tol = " << sci6(it.tol_val);
            }
            logging::info(true, os.str());
        }
    };

    logging::info(true, "");
    logging::info(true, "Post-checks");
    logging::up();

    for (const auto& it : items) {
        if (!it.show) continue;
        print_block(it);
    }

    logging::down();
}

}} // namespace fem::constraint
