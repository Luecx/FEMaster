/**
 * @file constraint_transformer_postcheck.cpp
 * @brief Implements diagnostic post-check utilities for `ConstraintTransformer`.
 *
 * @see src/constraints/transformer/constraint_transformer.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "constraint_transformer.h"

#include "../../core/logging.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace fem {
namespace constraint {

void ConstraintTransformer::post_check_static(const SparseMatrix& K,
                                              const DynamicVector& f,
                                              const DynamicVector& solution,
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

    const DynamicVector u    = recover_displacement(solution);
    const DynamicVector Cu_d = system_.C * u - system_.d;
    const Precision abs_Cu_d = Cu_d.norm();
    const Precision den_Cu_d = std::max<Precision>(1, system_.d.size() ? system_.d.norm() : 0);
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

        const DynamicVector residual           = K * u - f;
        const DynamicVector multipliers        = extract_lagrange_multipliers(solution);
        const DynamicVector scaled_multipliers = scale_lagrange_rows(multipliers);
        const DynamicVector kkt_u              = residual + system_.C.transpose() * scaled_multipliers;

        const Precision abs_kkt_u = kkt_u.norm();
        const Precision den_kkt_u = std::max<Precision>(1, std::max((K * u).norm(), f.norm()));
        const Precision rel_kkt_u = abs_kkt_u / den_kkt_u;

        items.push_back(Item{"lagrange equilibrium",
                             abs_kkt_u,
                             rel_kkt_u,
                             tol_reduced_rel,
                             rel_kkt_u <= tol_reduced_rel,
                             true});

        if (system_.equations > 0) {
            const Precision regularization = lagrange_regularization(K);
            DynamicVector kkt_c = scale_lagrange_rows(system_.C * u - system_.d);
            if (regularization > Precision(0)) {
                kkt_c.noalias() -= regularization * multipliers;
            }

            const Precision abs_kkt_c = kkt_c.norm();
            const Precision den_kkt_c = std::max<Precision>(
                1,
                std::max(
                    scale_lagrange_rows(system_.C * u).norm(),
                    scale_lagrange_rows(system_.d).norm()
                )
            );
            const Precision rel_kkt_c = abs_kkt_c / den_kkt_c;

            items.push_back(Item{"lagrange constraints",
                                 abs_kkt_c,
                                 rel_kkt_c,
                                 tol_constraint_rel,
                                 rel_kkt_c <= tol_constraint_rel,
                                 true});
        }

        print_items(items);
        return;
    }

    const DynamicVector resid = K * u - f;
    const Precision abs_resid = resid.norm();
    const Precision den_full = std::max<Precision>(1, f.size() ? f.norm() : 0);
    const Precision rel_resid = abs_resid / den_full;

    const Index nm = independent_dofs();
    DynamicVector Ttf(nm);
    DynamicVector TtKu(nm);
    DynamicVector red(nm);
    project_vector(f, Ttf);
    project_vector(K * u, TtKu);
    project_vector(resid, red);

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

    DynamicVector q0 = DynamicVector::Zero(nm);
    DynamicVector u_p = recover_displacement(q0);
    DynamicVector cup = system_.C * u_p - system_.d;
    const Precision abs = cup.norm();
    const Precision den = std::max<Precision>(1, system_.d.size() ? system_.d.norm() : 0);
    const Precision rel = abs / den;
    items.push_back(Item{"affine consistency (u_p)",
                         abs,
                         rel,
                         tol_constraint_rel,
                         rel <= tol_constraint_rel,
                         true});

    if (nm > 0) {
        const int samples = 3;
        Precision abs_acc = 0;
        Precision rel_acc = 0;
        const Precision CnormF = std::max<Precision>(1, system_.C.norm());

        for (int s = 0; s < samples; ++s) {
            DynamicVector q = DynamicVector::Random(nm);
            DynamicVector Tq;
            Tq = recover_increment(q);
            const Precision un = std::max<Precision>(Precision(1e-16), Tq.norm());
            DynamicVector v = system_.C * Tq;
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
