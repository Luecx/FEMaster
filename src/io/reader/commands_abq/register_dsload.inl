/**
 * @file register_dsload.inl
 * @brief Registers supported step-dependent Abaqus *DSLOAD definitions.
 *
 * Uniform pressure (`P`) and general traction (`TRVEC`) are retained as logical
 * surface-load records across steps. `OP=MOD` updates or adds records while
 * `OP=NEW` clears the complete active DSLOAD category before replacement data
 * are read.
 *
 * Abaqus identifies a general traction by surface, load type and specified load
 * direction. The unnormalized input direction is therefore retained for history
 * matching; FEMaster normalization is deferred until materialization. Repeated
 * definitions of one identity in a step remain separate additive conditions.
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
#include <stdexcept>
#include <string>

#include "../parser_abq.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers supported Abaqus surface-load history records.
 *
 * Pressure identity is surface plus `P`. General-traction identity additionally
 * includes the supplied direction components and orientation name because these
 * determine the physical traction direction. One propagated definition can be
 * redefined with `OP=MOD`; repeated definitions in the current step are additive.
 *
 * An empty `OP=NEW` block removes all active DSLOADs. Imaginary harmonic loads
 * remain unsupported.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser retaining active surface-load definitions.
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
            if (!state.step_active || !parser.active_loadcase()) {
                throw std::runtime_error("DSLOAD must appear after a supported procedure inside STEP");
            }
            if (state.procedure == "EIGENFREQ") {
                throw std::runtime_error("DSLOAD is not supported in a FREQUENCY step");
            }
            if (keys.has("REAL") && keys.has("IMAGINARY")) {
                throw std::runtime_error("DSLOAD REAL and IMAGINARY are mutually exclusive");
            }
            if (keys.has("IMAGINARY")) {
                throw std::runtime_error("DSLOAD IMAGINARY is not supported by the real-load harmonic solver");
            }

            const std::string op = keys.raw("OP");
            if (!state.dsload_op.empty() && state.dsload_op != op) {
                throw std::runtime_error("All DSLOAD cards in one STEP must use the same OP value");
            }
            if (state.dsload_op.empty()) {
                state.dsload_op = op;
                if (op == "NEW") {
                    state.dsloads.clear();
                }
            }

            *amplitude   = keys.has("AMPLITUDE")   ? keys.raw("AMPLITUDE")   : std::string{};
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            *follower    = keys.raw("FOLLOWER");

            if (!amplitude->empty() && !parser.model()._data->amplitudes.has(*amplitude)) {
                throw std::runtime_error("DSLOAD references unknown amplitude '" + *amplitude + "'");
            }
            if (!orientation->empty() && !parser.model()._data->coordinate_systems.has(*orientation)) {
                throw std::runtime_error("DSLOAD references unknown orientation '" + *orientation + "'");
            }
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
                    if (!model._data->surface_sets.has(surface)) {
                        throw std::runtime_error("DSLOAD references unknown surface '" + surface + "'");
                    }

                    std::array<fem::Precision, 3> stored_direction{
                        fem::Precision(0), fem::Precision(0), fem::Precision(0)
                    };
                    if (type == "P") {
                        if (!std::isnan(direction[0]) || !std::isnan(direction[1]) || !std::isnan(direction[2])) {
                            throw std::runtime_error("DSLOAD P accepts no traction direction components");
                        }
                    } else if (type == "TRVEC") {
                        if (std::isnan(direction[0]) || std::isnan(direction[1]) || std::isnan(direction[2])) {
                            throw std::runtime_error("DSLOAD TRVEC requires three direction components");
                        }

                        const fem::Vec3 traction_direction{direction[0], direction[1], direction[2]};
                        const fem::Precision norm = traction_direction.norm();
                        if (!(norm > fem::Precision(0)) || !std::isfinite(norm)) {
                            throw std::runtime_error("DSLOAD TRVEC direction must be finite and nonzero");
                        }
                        stored_direction = direction;
                    } else {
                        throw std::runtime_error("DSLOAD supports only P and TRVEC load labels");
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
                    if (matches > 1) {
                        throw std::runtime_error(
                            "Multiple active DSLOADs use the same surface, type and direction; use OP=NEW to redefine them"
                        );
                    }

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
