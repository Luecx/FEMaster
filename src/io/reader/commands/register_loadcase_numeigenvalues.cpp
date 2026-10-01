/**
 * @file register_loadcase_numeigenvalues.cpp
 * @brief Registers the requested number of eigenpairs for modal analyses.
 *
 * `NUMEIGENVALUES` reads one positive mode count inside a `LOADCASE` scope and
 * applies it to either linear buckling or eigenfrequency extraction. These are
 * the two FEMaster analyses whose result cardinality is determined by a
 * generalized eigenvalue solve.
 *
 * Spectral assembly and solver selection remain within the concrete load case;
 * this command only validates and stores the requested count.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../parser.h"

#include "../../../core/logging.h"
#include "../../../loadcase/linear_buckling.h"
#include "../../../loadcase/linear_eigenfreq.h"

#include <array>
#include <cmath>

namespace fem::io::reader::commands {

void register_loadcase_numeigenvalues(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("NUMEIGENVALUES", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Set number of eigenvalues for buckling/eigenfrequency loadcases.");

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<int>().name("COUNT").desc("Number of eigenvalues")
                )
                .bind([&parser](int count) {
                    logging::error(count > 0,
                        "NUMEIGENVALUES requires a positive integer");

                    auto* base = parser.active_loadcase();
                    logging::error(base != nullptr,
                        "NUMEIGENVALUES must appear inside *LOADCASE");

                    if (auto* lc = dynamic_cast<loadcase::LinearBuckling*>(base)) {
                        lc->num_eigenvalues = count;
                        return;
                    }
                    if (auto* lc = dynamic_cast<loadcase::LinearEigenfrequency*>(base)) {
                        lc->num_eigenvalues = count;
                        lc->use_eigenvalue_range = false;
                        return;
                    }

                    logging::error(false,
                        "NUMEIGENVALUES not supported for loadcase type ", base->type_name());
                })
            )
        );
    });

    registry.command("EIGENVALUERANGE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Set an open eigenvalue interval for eigenfrequency loadcases; bounds are lambda, not Hz.");

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<Precision, 2>().name("RANGE").desc("Minimum and maximum eigenvalue")
                )
                .bind([&parser](const std::array<Precision, 2>& bounds) {
                    auto* base = parser.active_loadcase();
                    logging::error(base != nullptr,
                        "EIGENVALUERANGE must appear inside *LOADCASE");

                    auto* lc = dynamic_cast<loadcase::LinearEigenfrequency*>(base);
                    logging::error(lc != nullptr,
                        "EIGENVALUERANGE not supported for loadcase type ", base->type_name());
                    logging::error(std::isfinite(bounds[0]) && std::isfinite(bounds[1]) && bounds[0] < bounds[1],
                        "EIGENVALUERANGE requires finite bounds with min < max");

                    lc->min_eigenvalue       = bounds[0];
                    lc->max_eigenvalue       = bounds[1];
                    lc->use_eigenvalue_range = true;
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
