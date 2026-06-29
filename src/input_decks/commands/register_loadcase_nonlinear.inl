// register_loadcase_nonlinear.inl — registers NONLINEAR controls within *LOADCASE

#include <stdexcept>

#include "../parser.h"

#include "../../loadcase/nonlinear_static.h"

namespace fem::input_decks::commands {

inline void register_loadcase_nonlinear(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("NONLINEAR", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Configure NONLINEARSTATIC increment and iteration controls.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("INCREMENTS").optional("10")
                .key("ADAPTIVE").optional("ON")
                .key("MAXITER").optional("20")
                .key("TOL").optional("1e-8")
                .key("REGULARIZE_ZERO_ROWS").optional("ON")
                .key("REGULARIZATION_ALPHA").optional("1e-4")
        );

        command.on_enter([&parser](const fem::dsl::Keys& keys) {
            auto* base = parser.active_loadcase();
            if (!base) {
                throw std::runtime_error("NONLINEAR must appear inside *LOADCASE");
            }

            auto* lc = dynamic_cast<loadcase::NonlinearStatic*>(base);
            if (!lc) {
                throw std::runtime_error("NONLINEAR is only supported for NONLINEARSTATIC loadcases");
            }

            lc->num_increments = keys.get<int>("INCREMENTS");
            lc->adaptive_increments = keys.get<bool>("ADAPTIVE");
            lc->max_iterations = keys.get<int>("MAXITER");
            lc->tolerance = keys.get<Precision>("TOL");
            lc->regularize_zero_stiffness_rows = keys.get<bool>("REGULARIZE_ZERO_ROWS");
            lc->zero_stiffness_regularization_alpha = keys.get<Precision>("REGULARIZATION_ALPHA");
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands
