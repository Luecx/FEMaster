// register_loadcase_numeigenvalues.inl — registers NUMEIGENVALUES for loadcases

#include <stdexcept>

#include "../parser.h"

#include "../../loadcase/linear_buckling.h"
#include "../../loadcase/linear_eigenfreq.h"

namespace fem::input_decks::commands {

inline void register_loadcase_numeigenvalues(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("NUMEIGENVALUES", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Set number of eigenvalues for buckling/eigenfrequency loadcases.");

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<int>().name("COUNT").desc("Number of eigenvalues")
                )
                .bind([&parser](int count) {
                    if (count <= 0) {
                        throw std::runtime_error("NUMEIGENVALUES requires a positive integer");
                    }

                    auto* base = parser.active_loadcase();
                    if (!base) {
                        throw std::runtime_error("NUMEIGENVALUES must appear inside *LOADCASE");
                    }

                    if (auto* lc = dynamic_cast<loadcase::LinearBuckling*>(base)) {
                        lc->num_eigenvalues = count;
                        return;
                    }
                    if (auto* lc = dynamic_cast<loadcase::LinearEigenfrequency*>(base)) {
                        lc->num_eigenvalues = count;
                        return;
                    }

                    throw std::runtime_error("NUMEIGENVALUES not supported for loadcase type " + parser.active_loadcase_type());
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands

