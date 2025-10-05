/**
 * @file constraint_transformer.cpp
 * @brief Implements the `ConstraintTransformer` facade utilities.
 *
 * The transformer orchestrates constraint assembly, map construction, and
 * provides diagnostic checks for static equilibrium solutions.
 *
 * @see src/constraints/constraint_transformer.h
 * @see src/constraints/builder/builder.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "constraint_transformer.h"

#include "../core/logging.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace fem {
namespace constraint {

/**
 * @copydoc ConstraintTransformer::ConstraintTransformer
 */
ConstraintTransformer::ConstraintTransformer(const Equations& eqs,
                                             const SystemDofIds& dofs,
                                             Index n_dofs,
                                             const BuildOptions& opt) {
    set_.equations = eqs;
    set_.opt = opt.set;
    set_.assemble(dofs, n_dofs);
    std::tie(map_, report_) = ConstraintBuilder::build(set_, opt.builder);
}

/**
 * @copydoc ConstraintTransformer::post_check_static
 */
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
        bool show_tol;
        Precision tol_val;
        bool pass;
        bool show;
    };

    DynamicVector Cu_d = set().C * u - set().d;
    const Precision abs_Cu_d = Cu_d.norm();
    const Precision den_Cu_d = std::max<Precision>(1, set().d.size() ? set().d.norm() : 0);
    const Precision rel_Cu_d = abs_Cu_d / (den_Cu_d > 0 ? den_Cu_d : 1);
    const bool pass_constraints = (rel_Cu_d <= tol_constraint_rel);

    DynamicVector resid = K * u - f;
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
    const bool pass_reduced = (rel_red <= tol_reduced_rel);

    std::vector<Item> items;
    items.reserve(5);

    items.push_back(Item{"constraints",
                         abs_Cu_d,
                         rel_Cu_d,
                         true,
                         tol_constraint_rel,
                         pass_constraints,
                         true});

    const bool evaluate_full = std::isfinite(tol_full_rel);
    const bool pass_full = evaluate_full ? (rel_resid <= tol_full_rel) : true;
    items.push_back(Item{"full residual",
                         abs_resid,
                         rel_resid,
                         true,
                         tol_full_rel,
                         pass_full,
                         evaluate_full});

    items.push_back(Item{"reduced equilibrium",
                         abs_red,
                         rel_red,
                         true,
                         tol_reduced_rel,
                         pass_reduced,
                         true});

    if (nm >= 0) {
        DynamicVector q0 = DynamicVector::Zero(nm);
        DynamicVector u_p = recover_u(q0);
        DynamicVector cup = set().C * u_p - set().d;
        const Precision abs = cup.norm();
        const Precision den = std::max<Precision>(1, set().d.size() ? set().d.norm() : 0);
        const Precision rel = abs / den;
        const bool pass = (rel <= tol_constraint_rel);
        items.push_back(Item{"affine consistency (u_p)",
                             abs,
                             rel,
                             true,
                             tol_constraint_rel,
                             pass,
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
        const bool pass = (rel_avg <= tol_constraint_rel);
        items.push_back(Item{"null-space (C*T=0)",
                             abs_avg,
                             rel_avg,
                             true,
                             tol_constraint_rel,
                             pass,
                             true});
    }

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

    const auto print_block = [&](const Item& item) {
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
              << rel_tag << sci6(item.rel_val);
        if (item.show_tol) {
            line2 << "  <= tol = " << sci6(item.tol_val);
        }
        logging::info(true, line2.str());
    };

    logging::info(true, "");
    logging::info(true, "Post-checks");
    logging::up();
    for (const auto& item : items) {
        if (!item.show) {
            continue;
        }
        print_block(item);
    }
    logging::down();
}

} // namespace constraint
} // namespace fem
