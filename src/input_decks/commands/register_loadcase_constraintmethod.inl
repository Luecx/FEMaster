// register_loadcase_constraintmethod.inl — registers CONSTRAINTMETHOD within *LOADCASE (LinearStatic and LinearStaticTopo)

#include <stdexcept>
#include <string>

#include "../parser.h"

#include "../../loadcase/linear_static.h"

namespace fem::input_decks::commands {

inline void register_loadcase_constraintmethod(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("CONSTRAINTMETHOD", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc(
            "Select constraint backend for LINEARSTATIC/LINEARSTATICTOPO: NULLSPACE or LAGRANGE.\n"
            "\n"
            "Constraint | Backend   | DIRECT       | INDIRECT\n"
            "NULLSPACE  | CPU MKL   | Yes          | Yes\n"
            "NULLSPACE  | CPU Eigen | Yes          | Yes\n"
            "NULLSPACE  | GPU       | Yes          | Yes\n"
            "NULLSPACE  | GPU cuDSS | Yes          | Yes\n"
            "LAGRANGE   | CPU MKL   | Yes          | No\n"
            "LAGRANGE   | CPU Eigen | Limited      | No\n"
            "LAGRANGE   | GPU       | No           | No\n"
            "LAGRANGE   | GPU cuDSS | Yes          | No"
        );

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("TYPE").required().allowed({"NULLSPACE", "LAGRANGE"})
        );

        command.on_enter([&parser](const fem::dsl::Keys& keys) {
            auto* base = parser.active_loadcase();
            if (!base) {
                throw std::runtime_error("CONSTRAINTMETHOD must appear inside *LOADCASE");
            }

            auto* lc = dynamic_cast<loadcase::LinearStatic*>(base);
            if (!lc) {
                throw std::runtime_error("CONSTRAINTMETHOD is only supported for LINEARSTATIC and LINEARSTATICTOPO loadcases");
            }

            const std::string type = keys.raw("TYPE");
            lc->constraint_method = (type == "LAGRANGE")
                ? constraint::ConstraintTransformer::Method::Lagrange
                : constraint::ConstraintTransformer::Method::NullSpace;
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands
