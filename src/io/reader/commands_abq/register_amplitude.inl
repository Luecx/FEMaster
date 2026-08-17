/**
 * @file register_amplitude.inl
 * @brief Registers tabular Abaqus *AMPLITUDE definitions.
 *
 * The initial Abaqus reader supports relative, step-time tabular amplitudes and
 * maps them directly to FEMaster's linearly interpolated `bc::Amplitude` model.
 * Abaqus permits up to four time/value pairs on one data line; the registration
 * expands every supplied pair into the same named FEMaster amplitude.
 *
 * Total-time, absolute-value and analytical/user-defined amplitude forms remain
 * unsupported because they require semantics beyond the existing FEMaster
 * scalar time-history representation.
 *
 * @see bc::Amplitude
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../bc/amplitude.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers Abaqus `DEFINITION=TABULAR` amplitude input.
 *
 * The supported form uses step time and relative values, matching the time
 * coordinate consumed by FEMaster transient load cases. Every line may contain
 * one to four `(time,value)` pairs. Incomplete pairs are rejected rather than
 * silently manufacturing a value.
 *
 * @param registry Stage-local DSL registry.
 * @param model FEMaster model receiving the named amplitude.
 */
inline void register_amplitude(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("AMPLITUDE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Define a relative step-time Abaqus tabular amplitude.");

        auto name = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME")
                    .required()
                    .doc("Amplitude identifier")
                .key("DEFINITION")
                    .optional("TABULAR")
                    .allowed({"TABULAR"})
                    .doc("Only TABULAR amplitudes are currently supported")
                .key("TIME")
                    .optional("STEPTIME")
                    .allowed({"STEPTIME"})
                    .doc("Amplitude time basis")
                .key("VALUE")
                    .optional("RELATIVE")
                    .allowed({"RELATIVE"})
                    .doc("Amplitude value interpretation")
        );

        command.on_enter([name, &model](const fem::io::dsl::Keys& keys) {
            *name = keys.raw("NAME");
            model.define_amplitude(*name, bc::Interpolation::Linear);
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 8>().name("DATA")
                        .desc("Up to four Abaqus time/value pairs")
                        .on_missing(std::numeric_limits<fem::Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<fem::Precision>::quiet_NaN())
                )
                .bind([&model, name](const std::array<fem::Precision, 8>& data) {
                    bool added = false;
                    for (std::size_t i = 0; i < data.size(); i += 2) {
                        const bool has_time  = !std::isnan(data[i]);
                        const bool has_value = !std::isnan(data[i + 1]);

                        if (!has_time && !has_value) {
                            continue;
                        }
                        if (has_time != has_value) {
                            throw std::runtime_error("AMPLITUDE requires complete time/value pairs");
                        }

                        model.add_amplitude_sample(*name, data[i], data[i + 1]);
                        added = true;
                    }

                    if (!added) {
                        throw std::runtime_error("AMPLITUDE data line contains no time/value pair");
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
