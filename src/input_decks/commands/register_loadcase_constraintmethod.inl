// register_loadcase_constraintmethod.inl — registers CONSTRAINTMETHOD within *LOADCASE (LinearStatic and LinearStaticTopo)

#include <stdexcept>
#include <string>

#include "../parser.h"

#include "../../loadcase/linear_static.h"
#include "../../loadcase/nonlinear_static.h"

namespace fem::input_decks::commands {

inline void register_loadcase_constraintmethod(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("CONSTRAINTMETHOD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc(
            "Select constraint backend for LINEARSTATIC/LINEARSTATICTOPO: NULLSPACE, LAGRANGE or ELIMINATION.\n"
            "\n"
            "Constraint | Backend   | DIRECT       | INDIRECT\n"
            "NULLSPACE  | CPU MKL   | Yes          | Yes\n"
            "NULLSPACE  | CPU Eigen | Yes          | Yes\n"
            "NULLSPACE  | GPU       | Yes          | Yes\n"
            "NULLSPACE  | GPU cuDSS | Yes          | Yes\n"
            "LAGRANGE   | CPU MKL   | Yes          | No\n"
            "LAGRANGE   | CPU Eigen | Limited      | No\n"
            "LAGRANGE   | GPU       | No           | No\n"
            "LAGRANGE   | GPU cuDSS | Yes          | No\n"
            "ELIMINATION| CPU MKL   | Yes          | Yes\n"
            "ELIMINATION| CPU Eigen | Yes          | Yes\n"
            "ELIMINATION| GPU       | Yes          | Yes\n"
            "ELIMINATION| GPU cuDSS | Yes          | Yes"
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

            auto* lc = dynamic_cast<loadcase::LinearStatic*>(base);
            auto* nlc = dynamic_cast<loadcase::NonlinearStatic*>(base);
            if (!lc && !nlc) {
                throw std::runtime_error("CONSTRAINTMETHOD is only supported for LINEARSTATIC, LINEARSTATICTOPO and NONLINEARSTATIC loadcases");
            }

            const std::string type = keys.raw("TYPE");
            auto method = constraint::ConstraintTransformer::Method::NullSpace;
            if (type == "LAGRANGE") {
                method = constraint::ConstraintTransformer::Method::Lagrange;
            } else if (type == "ELIMINATION") {
                method = constraint::ConstraintTransformer::Method::Elimination;
            }
            if (lc) {
                lc->constraint_method = method;
                return;
            }
            nlc->constraint_method = method;
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands
