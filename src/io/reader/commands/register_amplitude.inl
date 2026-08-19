/**
 * @file register_amplitude.inl
 * @brief Registers reusable scalar amplitudes.
 *
 * The root-level `AMPLITUDE` command builds named time-dependent scalar
 * functions from tabular samples. Supported interpolation modes are mapped to
 * `bc::Amplitude` and the completed object is stored in `ModelData` for later
 * use by transient, harmonic and other amplitude-aware loads.
 *
 * This file owns input validation and sample transfer only. Evaluation at a
 * physical analysis time remains the responsibility of the amplitude object.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../../../bc/amplitude.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_amplitude(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("AMPLITUDE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").required()
                .key("TYPE").optional("LINEAR").allowed({"STEP", "NEAREST", "LINEAR"})
        );
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            bc::Interpolation interpolation = bc::Interpolation::Linear;
            const std::string type = keys.raw("TYPE");
            if (type == "STEP") {
                interpolation = bc::Interpolation::Step;
            } else if (type == "NEAREST") {
                interpolation = bc::Interpolation::Nearest;
            }
            model.add_amplitude(std::make_shared<bc::Amplitude>(keys.raw("NAME"), interpolation));
        });
        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("TIME")
                    .one<fem::Precision>().name("VALUE")
                )
                .bind([&model](fem::Precision time, fem::Precision value) {
                    const auto amplitude = model._data->amplitudes.get();
                    logging::error(amplitude != nullptr,
                        "AMPLITUDE: no active amplitude is available");
                    amplitude->add_sample(time, value);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
