/**
 * @file register_dsload.inl
 * @brief Registers supported Abaqus *DSLOAD surface-load definitions.
 *
 * Uniform pressure (`P`) and general traction (`TRVEC`) are stored as logical
 * Abaqus surface-load records and propagated between analysis steps according to
 * `OP=MOD` and `OP=NEW`. The records are converted to FEMaster `PLoad` or `DLoad`
 * objects during `*END STEP` processing.
 *
 * General-traction identity contains the surface, load type, specified direction
 * and optional orientation. The input direction is retained until materialization
 * so Abaqus load-history matching uses the original definition.
 *
 * @see ParserAbqState
 * @see ParserAbqDSLoad
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>

#include "../parser_abq.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers the supported Abaqus surface-load syntax.
 *
 * Pressure records use `surface, P, magnitude`. General traction records use
 * `surface, TRVEC, magnitude, d1, d2, d3` and may reference an `ORIENTATION`.
 * `OP=MOD` preserves the active DSLOAD set and modifies a uniquely matching
 * propagated definition; `OP=NEW` clears the complete set before new records are
 * read. Repeated definitions with the same identity inside one step are additive.
 *
 * Imaginary harmonic load components are not supported. All DSLOAD cards within
 * one step must use the same `OP` value.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser containing the active surface-load state.
 */
inline void register_dsload(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("DSLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Define or modify active Abaqus pressure/general traction loads.");

        auto amplitude   = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();
        auto follower    = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("AMPLITUDE").optional().doc("Optional named Abaqus amplitude")
                .key("ORIENTATION").optional().doc("Optional traction coordinate system")
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
                    .doc("MOD preserves active DSLOADs; NEW replaces the complete DSLOAD category")
                .key("FOLLOWER").optional("YES").allowed({"YES", "NO"})
                    .doc("Abaqus general-traction follower setting")
                .flag("REAL").doc("Real harmonic load component")
                .flag("IMAGINARY").doc("Unsupported imaginary harmonic load component")
        );

        command.on_enter([&parser, amplitude, orientation, follower](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            logging::error(state.step_active && parser.active_loadcase(),
                "DSLOAD must appear after a supported procedure inside STEP");
            logging::error(state.procedure != "EIGENFREQ",
                "DSLOAD is not supported in a FREQUENCY step");
            logging::error(!(keys.has("REAL") && keys.has("IMAGINARY")),
                "DSLOAD REAL and IMAGINARY are mutually exclusive");
            logging::error(!keys.has("IMAGINARY"),
                "DSLOAD IMAGINARY is not supported by the real-load harmonic solver");

            const std::string op = keys.raw("OP");
            logging::error(state.dsload_op.empty() || state.dsload_op == op,
                "All DSLOAD cards in one STEP must use the same OP value");
            if (state.dsload_op.empty()) {
                state.dsload_op = op;
                if (op == "NEW") {
                    state.dsloads.clear();
                }
            }

            *amplitude   = keys.has("AMPLITUDE")   ? keys.raw("AMPLITUDE")   : std::string{};
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            *follower    = keys.raw("FOLLOWER");

            logging::error(amplitude->empty() || parser.model()._data->amplitudes.has(*amplitude),
                "DSLOAD references unknown amplitude '", *amplitude, "'");
            logging::error(orientation->empty() || parser.model()._data->coordinate_systems.has(*orientation),
                "DSLOAD references unknown orientation '", *orientation, "'");
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("SURFACE").desc("Abaqus surface name")
                    .one<std::string>().name("TYPE").desc("Supported labels P or TRVEC")
                    .one<fem::Precision>().name("MAGNITUDE").desc("Reference load magnitude")
                    .fixed<fem::Precision, 3>().name("DIRECTION").desc("TRVEC direction components")
                        .on_missing(std::numeric_limits<fem::Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<fem::Precision>::quiet_NaN())
                )
                .bind([&parser, amplitude, orientation, follower](const std::string& surface,
                                                                  const std::string& type,
                                                                  fem::Precision magnitude,
                                                                  const std::array<fem::Precision, 3>& direction) {
                    auto& model = parser.model();
                    auto& state = parser.abaqus_state();
                    logging::error(model._data->surface_sets.has(surface),
                        "DSLOAD references unknown surface '", surface, "'");

                    std::array<fem::Precision, 3> stored_direction{
                        fem::Precision(0), fem::Precision(0), fem::Precision(0)
                    };
                    if (type == "P") {
                        logging::error(std::isnan(direction[0])
                                    && std::isnan(direction[1])
                                    && std::isnan(direction[2]),
                            "DSLOAD P accepts no traction direction components");
                    } else if (type == "TRVEC") {
                        logging::error(!std::isnan(direction[0])
                                    && !std::isnan(direction[1])
                                    && !std::isnan(direction[2]),
                            "DSLOAD TRVEC requires three direction components");

                        const fem::Vec3 traction_direction{direction[0], direction[1], direction[2]};
                        const fem::Precision norm = traction_direction.norm();
                        logging::error(norm > fem::Precision(0) && std::isfinite(norm),
                            "DSLOAD TRVEC direction must be finite and nonzero");
                        stored_direction = direction;
                    } else {
                        logging::error(false,
                            "DSLOAD supports only P and TRVEC load labels");
                    }

                    const auto same_identity = [&](const ParserAbqDSLoad& item) {
                        if (item.surface != surface || item.type != type) {
                            return false;
                        }
                        if (type == "P") {
                            return true;
                        }
                        return item.direction == stored_direction && item.orientation == *orientation;
                    };

                    auto first = state.dsloads.end();
                    int matches = 0;
                    bool defined_this_step = false;

                    for (auto it = state.dsloads.begin(); it != state.dsloads.end(); ++it) {
                        if (!same_identity(*it)) {
                            continue;
                        }
                        if (first == state.dsloads.end()) {
                            first = it;
                        }
                        ++matches;
                        defined_this_step = defined_this_step || it->modified_step == state.step_index;
                    }

                    if (defined_this_step || matches == 0) {
                        state.dsloads.push_back(ParserAbqDSLoad{
                            surface,
                            type,
                            magnitude,
                            fem::Precision(0),
                            stored_direction,
                            *amplitude,
                            *orientation,
                            *follower,
                            state.step_index
                        });
                        return;
                    }

                    logging::error(matches == 1,
                        "Multiple active DSLOADs use the same surface, type and direction; use OP=NEW to redefine them");

                    first->previous_magnitude = first->magnitude;
                    first->magnitude          = magnitude;
                    first->amplitude          = *amplitude;
                    first->follower           = *follower;
                    first->modified_step      = state.step_index;
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
