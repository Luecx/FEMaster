/**
 * @file register_dload.inl
 * @brief Registers supported step-dependent Abaqus *DLOAD definitions.
 *
 * The current element-based distributed-load subset supports `GRAV`. Definitions
 * are retained logically across steps so `OP=MOD` can update an existing target
 * and `OP=NEW` can replace the complete active DLOAD category without carrying
 * any mechanical solver state between FEMaster load cases.
 *
 * @see ParserAbqState
 * @see ParserAbqDLoad
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <algorithm>
#include <array>
#include <memory>
#include <stdexcept>
#include <string>

#include "../parser_abq.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers Abaqus gravity-load history records.
 *
 * The logical identity consists of the original element/element-set target and
 * the Abaqus load type. Materialization to FEMaster `VLoad` objects is deferred
 * until `*END STEP`. `OP` must be consistent across all `*DLOAD` cards in one
 * step.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser retaining active distributed-load definitions.
 */
inline void register_dload(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("DLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Define or modify active Abaqus GRAV distributed loads.");

        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("AMPLITUDE").optional().doc("Optional named Abaqus amplitude")
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
                    .doc("MOD preserves active DLOADs; NEW replaces the complete DLOAD category")
                .flag("REAL").doc("Real harmonic load component")
                .flag("IMAGINARY").doc("Unsupported imaginary harmonic load component")
        );

        command.on_enter([&parser, amplitude](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || !parser.active_loadcase()) {
                throw std::runtime_error("DLOAD must appear after a supported procedure inside STEP");
            }
            if (state.procedure == "EIGENFREQ") {
                throw std::runtime_error("DLOAD is not supported in a FREQUENCY step");
            }
            if (keys.has("REAL") && keys.has("IMAGINARY")) {
                throw std::runtime_error("DLOAD REAL and IMAGINARY are mutually exclusive");
            }
            if (keys.has("IMAGINARY")) {
                throw std::runtime_error("DLOAD IMAGINARY is not supported by the real-load harmonic solver");
            }

            const std::string op = keys.raw("OP");
            if (!state.dload_op.empty() && state.dload_op != op) {
                throw std::runtime_error("All DLOAD cards in one STEP must use the same OP value");
            }
            if (state.dload_op.empty()) {
                state.dload_op = op;
                if (op == "NEW") {
                    state.dloads.clear();
                }
            }

            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            if (!amplitude->empty() && !parser.model()._data->amplitudes.has(*amplitude)) {
                throw std::runtime_error("DLOAD references unknown amplitude '" + *amplitude + "'");
            }
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Element set, element id or blank for all elements")
                        .on_empty(std::string{"EALL"})
                    .one<std::string>().name("TYPE").desc("Supported load label GRAV")
                    .one<fem::Precision>().name("MAGNITUDE").desc("Gravity magnitude")
                    .fixed<fem::Precision, 3>().name("DIRECTION").desc("Gravity vector components")
                )
                .bind([&parser, amplitude](const std::string& target,
                                           const std::string& type,
                                           fem::Precision magnitude,
                                           const std::array<fem::Precision, 3>& direction) {
                    if (type != "GRAV") {
                        throw std::runtime_error(
                            "DLOAD currently supports only GRAV; use SURFACE + DSLOAD for surface pressure/traction"
                        );
                    }

                    auto& state = parser.abaqus_state();
                    auto record = std::find_if(
                        state.dloads.begin(), state.dloads.end(),
                        [&](const ParserAbqDLoad& item) {
                            return item.target == target && item.type == type;
                        }
                    );

                    const ParserAbqDLoad value{
                        target,
                        type,
                        magnitude,
                        direction,
                        *amplitude,
                        state.step_index
                    };

                    if (record == state.dloads.end()) {
                        state.dloads.push_back(value);
                    } else {
                        *record = value;
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
