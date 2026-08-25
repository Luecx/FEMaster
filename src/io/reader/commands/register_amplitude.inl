/**
 * @file register_amplitude.inl
 * @brief Registers the AMPLITUDE input command.
 *
 * The root-level `AMPLITUDE` command creates named time-dependent scalar
 * functions from tabular `(time, value)` samples. The interpolation mode is
 * mapped to `bc::Interpolation` and the resulting `bc::Amplitude` is stored
 * in the model for later use by amplitude-aware loads and analyses.
 *
 * The active amplitude pointer is maintained by the reader and identifies the
 * object currently receiving tabular samples. Evaluation at a physical
 * analysis time remains the responsibility of `bc::Amplitude`.
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

namespace dsl = fem::io::dsl;

/**
 * @brief Registers the root-level `AMPLITUDE` command.
 *
 * The command creates the named amplitude when its keyword line is entered.
 * Subsequent data lines append `(time, value)` samples to the active amplitude
 * maintained by the reader.
 */
inline void register_amplitude(dsl::Registry&      registry,
                               model::Model&       model) {

    bc::Amplitude::Ptr amplitude = std::make_shared<bc::Amplitude>(nullptr);

    registry.command("AMPLITUDE", [&](dsl::Command& command) {
        // Restrict AMPLITUDE to the root scope
        command.allow_if(dsl::Condition::parent_is("ROOT"));

        // Define the amplitude name and interpolation mode
        command.keyword(
            dsl::KeywordSpec::make()
                .key("NAME").required()
                .key("TYPE").optional("LINEAR").allowed({"LINEAR", "STEP", "NEAREST"})
        );

        // Create the amplitude before processing its tabular samples
        command.on_enter([&model, &amplitude](const dsl::Keys& keys) {
            bc::Interpolation interpolation = bc::Interpolation::Linear;

            const std::string type = keys.raw("TYPE");

            // Map the textual interpolation mode to the amplitude representation
            if (type == "STEP") {
                interpolation = bc::Interpolation::Step;
            } else if (type == "NEAREST") {
                interpolation = bc::Interpolation::Nearest;
            }

            // Keep the new amplitude active while its data lines are processed
            amplitude = std::make_shared<bc::Amplitude>(keys.raw("NAME"), interpolation);

            // Store the amplitude permanently in the model
            model.add_amplitude(amplitude);
        });

        // Read one (time, value) sample from every following data line
        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .one<fem::Precision>().name("TIME")
                    .one<fem::Precision>().name("VALUE")
                )
                .bind([&amplitude](fem::Precision time, fem::Precision value) {
                    // Ensure that the keyword created an amplitude before data is consumed
                    logging::error(amplitude != nullptr,
                        "AMPLITUDE: no active amplitude is available");

                    // Append the tabular sample to the active amplitude
                    amplitude->add_sample(time, value);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands