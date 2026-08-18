// register_loadcase_numeigenvalues.inl — registers NUMEIGENVALUES for loadcases

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
                        "NUMEIGENVALUES not supported for loadcase type ", parser.active_loadcase_type());
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands

