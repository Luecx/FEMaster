/**
 * @file register_amplitude.inl
 * @brief Registers Abaqus tabular amplitudes.
 *
 * The root-level Abaqus `AMPLITUDE` command converts pairs of independent and
 * dependent values into a named FEMaster `bc::Amplitude`. Supported Abaqus
 * interpolation and time-basis options are validated before the completed
 * function is stored in `ModelData`.
 *
 * Later step-load commands resolve the name and either retain the amplitude for
 * dynamic evaluation or sample it according to the active static procedure.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <cmath>
#include <limits>

#include "../../../bc/amplitude.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

inline void register_amplitude(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("AMPLITUDE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").required()
                .key("DEFINITION").optional("TABULAR").allowed({"TABULAR"})
                .key("TIME").optional("STEPTIME").allowed({"STEPTIME"})
                .key("VALUE").optional("RELATIVE").allowed({"RELATIVE"})
        );
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            model.add_amplitude(std::make_shared<bc::Amplitude>(
                keys.raw("NAME"), bc::Interpolation::Linear));
        });
        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<Precision, 8>().name("DATA")
                        .on_missing(std::numeric_limits<Precision>::quiet_NaN())
                        .on_empty(std::numeric_limits<Precision>::quiet_NaN())
                )
                .bind([&model](const std::array<Precision, 8>& data) {
                    const auto amplitude = model._data->amplitudes.get();
                    logging::error(amplitude != nullptr,
                        "AMPLITUDE: no active amplitude is available");

                    bool added = false;
                    for (std::size_t i = 0; i < data.size(); i += 2) {
                        const bool has_time  = !std::isnan(data[i]);
                        const bool has_value = !std::isnan(data[i + 1]);
                        if (!has_time && !has_value) continue;
                        logging::error(has_time && has_value,
                            "AMPLITUDE: incomplete time/value pair");
                        amplitude->add_sample(data[i], data[i + 1]);
                        added = true;
                    }
                    logging::error(added,
                        "AMPLITUDE: data line contains no time/value pair");
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
