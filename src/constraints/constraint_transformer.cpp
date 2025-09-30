/******************************************************************************
 * @file constraint_transformer.cpp
 * @brief Diagnostics & checklist for static post-checks on ConstraintTransformer.
 ******************************************************************************/

#include "constraint_transformer.h"
#include "../core/logging.h"

#include <cmath>       // std::isfinite
#include <iomanip>     // std::setw, std::setprecision, std::scientific
#include <sstream>     // std::ostringstream
#include <vector>
#include <string>
#include <algorithm>
#include <limits>

namespace fem { namespace constraint {

void ConstraintTransformer::post_check_static(const SparseMatrix& K,
                                              const DynamicVector& f,
                                              const DynamicVector& u,
                                              Precision tol_constraint_rel,
                                              Precision tol_reduced_rel,
                                              Precision tol_full_rel) const
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

    // -------------------- compute core quantities --------------------
    // constraints @ u
    DynamicVector Cu_d = set().C * u - set().d;
    const Precision abs_Cu_d   = Cu_d.norm();
    const Precision den_Cu_d   = std::max<Precision>(1, set().d.size() ? set().d.norm() : 0);
    const Precision rel_Cu_d   = abs_Cu_d / (den_Cu_d > 0 ? den_Cu_d : 1);
    const bool pass_constraints = (rel_Cu_d <= tol_constraint_rel);

    // full residual
    DynamicVector resid = K * u - f;
    const Precision abs_resid = resid.norm();
    const Precision den_full  = std::max<Precision>(1, f.size() ? f.norm() : 0);
    const Precision rel_resid = abs_resid / den_full;

    // reduced residual
    const Index nm = n_master();
    DynamicVector Ttf(nm), TtKu(nm), red(nm);
    apply_Tt(f,      Ttf);
    apply_Tt(K * u,  TtKu);
    apply_Tt(resid,  red);

    const Precision abs_red   = red.norm();
    const Precision den_red   = std::max<Precision>(1, std::max(Ttf.norm(), TtKu.norm()));
    const Precision rel_red   = abs_red / den_red;
    const bool pass_reduced   = (rel_red <= tol_reduced_rel);

    // -------------------- build checklist items ----------------------
    std::vector<Item> items;

    items.push_back(Item{
        "constraints",
        abs_Cu_d, rel_Cu_d,
        /*show_tol*/true, tol_constraint_rel,
        pass_constraints,
        /*show*/true
    });

    items.push_back(Item{
        "full residual",
        abs_resid, rel_resid,
        /*show_tol*/true, tol_full_rel,
        std::isfinite(tol_full_rel) ? (rel_resid <= tol_full_rel) : true,
        /*show*/std::isfinite(tol_full_rel)
    });

    items.push_back(Item{
        "reduced equilibrium",
        abs_red, rel_red,
        /*show_tol*/true, tol_reduced_rel,
        pass_reduced,
        /*show*/true
    });

    // Extra check: affine consistency (C * u_p ≈ d)
    {
        if (nm >= 0) {
            DynamicVector q0 = DynamicVector::Zero(nm);
            DynamicVector u_p = recover_u(q0);                 // u_p + T*0 = u_p
            DynamicVector cup = set().C * u_p - set().d;

            const Precision abs  = cup.norm();
            const Precision den  = std::max<Precision>(1, set().d.size() ? set().d.norm() : 0);
            const Precision rel  = abs / den;
            const bool      pass = (rel <= tol_constraint_rel);

            items.push_back(Item{
                "affine consistency (u_p)",
                abs, rel,
                /*show_tol*/true, tol_constraint_rel,
                pass,
                /*show*/true
            });
        }
    }

    // Extra check: null-space property (C*T ≈ 0) via 3 random samples
    {
        if (nm > 0) {
            const int S = 3;
            Precision abs_acc = 0;
            Precision rel_acc = 0;

            const Precision CnormF = std::max<Precision>(1, set().C.norm());

            for (int s = 0; s < S; ++s) {
                DynamicVector q  = DynamicVector::Random(nm);
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
            const bool      pass    = (rel_avg <= tol_constraint_rel);

            items.push_back(Item{
                "null-space (C*T=0)",
                abs_avg, rel_avg,
                /*show_tol*/true, tol_constraint_rel,
                pass,
                /*show*/true
            });
        }
    }

    // -------------------- dynamic label width + print ----------------
    int label_w = 22;
    for (const auto& it : items)
        if (it.show) label_w = std::max<int>(label_w, (int)it.label.size());

    const std::string bullet    = "  ";
    const std::string sp        = " ";
    const std::string colon_sep = " : ";
    const std::string abs_tag   = "abs = ";
    const std::string rel_tag   = "rel = ";
    const std::string tag_pad(6, ' '); // width of "[PASS]" or "[FAIL]"

    auto print_block = [&](const Item& it)
    {
        // Line 1
        {
            std::ostringstream os;
            os << bullet
               << (it.pass ? "[PASS]" : "[FAIL]") << sp
               << std::left << std::setw(label_w) << it.label
               << colon_sep
               << abs_tag << sci6(it.abs_val);
            logging::info(true, os.str());
        }
        // Line 2
        {
            std::ostringstream os;
            os << bullet
               << tag_pad << sp
               << std::left << std::setw(label_w) << ""
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
