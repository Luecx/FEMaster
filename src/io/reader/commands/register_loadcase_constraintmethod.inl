/**
 * @file register_loadcase_constraintmethod.inl
 * @brief Registers constraint-transformation selection for structural analyses.
 *
 * `CONSTRAINTMETHOD` selects null-space projection, Lagrange multipliers or
 * elimination for the active supported load case. The command translates the
 * deck token into `ConstraintTransformer::Method` and applies it only to
 * analyses that expose a compatible constraint backend.
 *
 * Solver/backend compatibility remains documented by the command grammar and
 * is enforced by the corresponding load-case implementation during execution.
 * Unsupported active load-case types produce a diagnostic containing their
 * canonical `type_name()`.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

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
            "LINEARHARMONIC currently accepts only NULLSPACE during execution.\n"
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
            logging::error(base != nullptr,
                "CONSTRAINTMETHOD must appear inside *LOADCASE");

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

            logging::error(false,
                "CONSTRAINTMETHOD not supported for loadcase type " + base->type_name());
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
