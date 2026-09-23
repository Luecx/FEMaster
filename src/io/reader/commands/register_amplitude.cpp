/**
 * @file register_amplitude.cpp
 * @brief Registers shared tabular AMPLITUDE input for both readers.
 *
 * The root-level `AMPLITUDE` command creates named time-dependent scalar
 * functions from up to four `(time, value)` pairs per data line. FEMaster's
 * `TYPE` selects the interpolation mode; supported Abaqus `DEFINITION`, `TIME`
 * and `VALUE` options describe the same tabular function.
 *
 * The active amplitude pointer is maintained by the reader and identifies the
 * object currently receiving tabular samples. Evaluation at a physical
 * analysis time remains the responsibility of `bc::Amplitude`.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../../../bc/amplitude.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>

namespace fem::io::reader::commands {

namespace dsl = fem::io::dsl;

/**
 * @brief Registers the root-level `AMPLITUDE` command.
 *
 * The command creates the named amplitude when its keyword line is entered.
 * Subsequent data lines append one to four `(time, value)` samples to the active
 * amplitude. Both readers use the same grammar and `bc::Amplitude` state.
 * Omitted pairs are ignored, while incomplete pairs are rejected.
 *
 * @param registry Stage-local DSL registry.
 * @param model FEMaster model receiving the named amplitude.
 */
void register_amplitude(dsl::Registry& registry, model::Model& model) {
    auto amplitude = std::make_shared<bc::Amplitude::Ptr>();

    registry.command("AMPLITUDE", [&](dsl::Command& command) {
        // Restrict AMPLITUDE to the root scope
        command.allow_if(dsl::Condition::parent_is("ROOT"));

        // Accept FEMaster interpolation and supported Abaqus tabular options.
        command.keyword(
            dsl::KeywordSpec::make()
                .key("NAME").required()
                .key("TYPE").optional("LINEAR").allowed({"LINEAR", "STEP", "NEAREST"})
                .key("DEFINITION").optional("TABULAR").allowed({"TABULAR"})
                .key("TIME").optional("STEPTIME").allowed({"STEPTIME"})
                .key("VALUE").optional("RELATIVE").allowed({"RELATIVE"})
        );

        // Create the amplitude before processing its tabular samples
        command.on_enter([&model, amplitude](const dsl::Keys& keys) {
            bc::Interpolation interpolation = bc::Interpolation::Linear;

            const std::string type = keys.raw("TYPE");

            // Map the textual interpolation mode to the amplitude representation
            if (type == "STEP") {
                interpolation = bc::Interpolation::Step;
            } else if (type == "NEAREST") {
                interpolation = bc::Interpolation::Nearest;
            }

            // Keep the new amplitude active while its data lines are processed
            *amplitude = std::make_shared<bc::Amplitude>(keys.raw("NAME"), interpolation);

            // Store the amplitude permanently in the model
            model.add_amplitude(*amplitude);
        });

        // Decode each physical line as up to four independent time/value pairs.
        command.variant(dsl::Variant::make()
            .segment(dsl::Segment::make()
                .range(dsl::LineRange{}.min(1))
                .pattern(dsl::Pattern::make()
                    .fixed<Precision, 8>().name("DATA").desc("Up to four time/value pairs")
                        .on_missing(std::numeric_limits<Precision>::quiet_NaN())
                        .on_empty(std::numeric_limits<Precision>::quiet_NaN())
                )
                .bind([amplitude](const std::array<Precision, 8>& data) {
                    // Ensure that the keyword created an amplitude before data is consumed
                    logging::error(*amplitude != nullptr,
                        "AMPLITUDE: no active amplitude is available");

                    // Preserve pair order; absent pairs carry NaN in both positions.
                    bool added = false;
                    for (std::size_t i = 0; i < data.size(); i += 2) {
                        const bool has_time  = !std::isnan(data[i]);
                        const bool has_value = !std::isnan(data[i + 1]);
                        if (!has_time && !has_value) continue;
                        logging::error(has_time && has_value,
                            "AMPLITUDE: incomplete time/value pair");
                        (*amplitude)->add_sample(data[i], data[i + 1]);
                        added = true;
                    }
                    logging::error(added,
                        "AMPLITUDE: data line contains no time/value pair");
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
