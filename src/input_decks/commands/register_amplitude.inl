// register_amplitude.inl — DSL registration for *AMPLITUDE

#include <memory>
#include <string>

#include "../../bc/amplitude.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_amplitude(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("AMPLITUDE", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Define a reusable scalar time history. Each data line specifies a time-value pair. "
            "The TYPE keyword controls interpolation (STEP, NEAREST, LINEAR). Loads referencing the amplitude "
            "will scale their components by the interpolated value at the current analysis time."
        );

        auto name = std::make_shared<std::string>();
        auto interpolation = std::make_shared<bc::Interpolation>(bc::Interpolation::Linear);

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("NAME").required().doc("Amplitude identifier")
                .key("TYPE").optional("LINEAR")
                    .doc("Interpolation scheme: STEP, NEAREST, or LINEAR (default)")
                    .allowed({"STEP", "NEAREST", "LINEAR"})
        );

        command.on_enter([name, interpolation, &model](const fem::dsl::Keys& keys) {
            *name = keys.raw("NAME");
            const std::string type_token = keys.raw("TYPE");
            if (type_token == "STEP" || type_token == "step" || type_token == "Step") {
                *interpolation = bc::Interpolation::Step;
            } else if (type_token == "NEAREST" || type_token == "nearest" || type_token == "Nearest") {
                *interpolation = bc::Interpolation::Nearest;
            } else {
                *interpolation = bc::Interpolation::Linear;
            }
            model.define_amplitude(*name, *interpolation);
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<fem::Precision>().name("TIME").desc("Time coordinate")
                    .one<fem::Precision>().name("VALUE").desc("Amplitude value at TIME")
                )
                .bind([&model, name](fem::Precision time, fem::Precision value) {
                    model.add_amplitude_sample(*name, time, value);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
