// register_loadcase_constraintmethod.inl — registers CONSTRAINTMETHOD within *LOADCASE

#include <stdexcept>
#include <string>

#include "../parser.h"

#include "../../../loadcase/linear_harmonic.h"
#include "../../../loadcase/linear_static.h"
#include "../../../loadcase/nonlinear_static.h"

namespace fem::io::reader::commands {

inline void register_loadcase_constraintmethod(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("CONSTRAINTMETHOD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc(
            "Select constraint backend for supported structural loadcases: NULLSPACE, LAGRANGE or ELIMINATION.\n"
            "LINEARHARMONIC currently accepts only NULLSPACE during execution."
        );

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE").required().allowed({"NULLSPACE", "LAGRANGE", "ELIMINATION"})
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto* base = parser.active_loadcase();
            if (!base) {
                throw std::runtime_error("CONSTRAINTMETHOD must appear inside *LOADCASE");
            }

            const std::string type = keys.raw("TYPE");
            auto method = constraint::ConstraintTransformer::Method::NullSpace;
            if (type == "LAGRANGE") {
                method = constraint::ConstraintTransformer::Method::Lagrange;
            } else if (type == "ELIMINATION") {
                method = constraint::ConstraintTransformer::Method::Elimination;
            }

            if (auto* lc = dynamic_cast<loadcase::LinearStatic*>(base)) {
                lc->constraint_method = method;
                return;
            }
            if (auto* lc = dynamic_cast<loadcase::NonlinearStatic*>(base)) {
                lc->constraint_method = method;
                return;
            }
            if (auto* lc = dynamic_cast<loadcase::LinearHarmonic*>(base)) {
                lc->constraint_method = method;
                return;
            }

            throw std::runtime_error(
                "CONSTRAINTMETHOD not supported for loadcase type " +
                parser.active_loadcase_type());
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
