/**
 * @file register_loadcase_numeigenvalues.inl
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

#include "../parser.h"

#include "../../../core/logging.h"
#include "../../../loadcase/linear_buckling.h"
#include "../../../loadcase/linear_eigenfreq.h"

namespace fem::io::reader::commands {

inline void register_loadcase_numeigenvalues(fem::io::dsl::Registry& registry, Parser& parser) {
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
                        return;
                    }

                    logging::error(false,
                        "NUMEIGENVALUES not supported for loadcase type ", base->type_name());
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
