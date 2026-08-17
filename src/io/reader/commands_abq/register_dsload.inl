/**
 * @file register_dsload.inl
 * @brief Registers supported Abaqus *DSLOAD surface loads.
 *
 * Uniform surface pressure (`P`) maps to FEMaster `PLoad`, while uniform general
 * surface traction (`TRVEC`) maps to `DLoad`. Abaqus traction direction vectors
 * are normalized before multiplication by the supplied reference magnitude.
 * Optional `ORIENTATION` names reuse FEMaster coordinate systems created by
 * `*ORIENTATION` and are stored directly on the traction load.
 *
 * Named amplitudes follow the same supported procedure semantics as concentrated
 * loads: transient histories remain attached, linear static/buckling loads are
 * reduced to their final step value, and unsupported nonlinear/harmonic history
 * interpretations are rejected.
 *
 * @see bc::DLoad
 * @see bc::PLoad
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include "../parser_abq.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers uniform Abaqus surface pressure and general traction records.
 *
 * `P` accepts `surface, P, magnitude`. `TRVEC` additionally requires three
 * direction components and normalizes that direction as Abaqus does. Follower
 * traction semantics are intentionally rejected because FEMaster's vector
 * surface load currently stores a reference/global or named-coordinate-system
 * direction without an Abaqus follower flag.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser providing current step state.
 */
inline void register_dsload(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("DSLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Apply uniform Abaqus pressure or general traction to a named surface.");

        auto amplitude   = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("AMPLITUDE").optional().doc("Optional named Abaqus amplitude")
                .key("ORIENTATION").optional().doc("Optional traction coordinate system")
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
                    .doc("Accepted for the independent current-step collector")
                .flag("FOLLOWER").doc("Unsupported Abaqus follower traction flag")
        );

        command.on_enter([&parser, amplitude, orientation](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || !parser.active_loadcase()) {
                throw std::runtime_error("DSLOAD must appear after a supported procedure inside STEP");
            }
            if (keys.has("FOLLOWER")) {
                throw std::runtime_error("DSLOAD FOLLOWER is not supported");
            }

            *amplitude   = keys.has("AMPLITUDE")   ? keys.raw("AMPLITUDE")   : std::string{};
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};

            if (!orientation->empty() && !parser.model()._data->coordinate_systems.has(*orientation)) {
                throw std::runtime_error("DSLOAD references unknown orientation '" + *orientation + "'");
            }

            parser.model()._data->load_cols.activate(state.load_collector);
            state.load_collector_used = true;
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("SURFACE").desc("Abaqus surface name")
                    .one<std::string>().name("TYPE").desc("Supported labels P or TRVEC")
                    .one<fem::Precision>().name("MAGNITUDE").desc("Reference load magnitude")
                    .fixed<fem::Precision, 3>().name("DIRECTION").desc("TRVEC direction components")
                        .on_missing(std::numeric_limits<fem::Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<fem::Precision>::quiet_NaN())
                )
                .bind([&parser, amplitude, orientation](const std::string& surface,
                                                        const std::string& type,
                                                        fem::Precision magnitude,
                                                        const std::array<fem::Precision, 3>& direction) {
                    auto& model = parser.model();
                    auto& state = parser.abaqus_state();
                    if (!model._data->surface_sets.has(surface)) {
                        throw std::runtime_error("DSLOAD references unknown surface '" + surface + "'");
                    }

                    std::string stored_amplitude;
                    fem::Precision scale = fem::Precision(1);
                    if (!amplitude->empty()) {
                        if (!model._data->amplitudes.has(*amplitude)) {
                            throw std::runtime_error("DSLOAD references unknown amplitude '" + *amplitude + "'");
                        }

                        if (state.procedure == "LINEARTRANSIENT") {
                            stored_amplitude = *amplitude;
                        } else if (state.procedure == "LINEARSTATIC" ||
                                   state.procedure == "LINEARBUCKLING") {
                            scale = model._data->amplitudes.get(*amplitude)->evaluate(state.step_period);
                        } else if (state.procedure == "NONLINEARSTATIC" ||
                                   state.procedure == "STATIC_RIKS") {
                            throw std::runtime_error(
                                "DSLOAD AMPLITUDE is not supported for nonlinear static/Riks proportional loading"
                            );
                        } else if (state.procedure == "LINEARHARMONIC") {
                            throw std::runtime_error(
                                "DSLOAD AMPLITUDE is not supported for harmonic response"
                            );
                        }
                    } else if (state.procedure == "LINEARTRANSIENT" &&
                               !state.step_amplitude.empty() &&
                               state.step_amplitude != "STEP") {
                        stored_amplitude = "__ABQ_STEP_" + std::to_string(state.step_index) + "_DEFAULT_AMPLITUDE";
                    }

                    if (type == "P") {
                        if (!std::isnan(direction[0]) || !std::isnan(direction[1]) || !std::isnan(direction[2])) {
                            throw std::runtime_error("DSLOAD P accepts no traction direction components");
                        }
                        model.add_pload(surface, scale * magnitude, stored_amplitude);
                        return;
                    }

                    if (type != "TRVEC") {
                        throw std::runtime_error("DSLOAD supports only P and TRVEC load labels");
                    }
                    if (std::isnan(direction[0]) || std::isnan(direction[1]) || std::isnan(direction[2])) {
                        throw std::runtime_error("DSLOAD TRVEC requires three direction components");
                    }

                    fem::Vec3 unit_direction{direction[0], direction[1], direction[2]};
                    const fem::Precision norm = unit_direction.norm();
                    if (!(norm > fem::Precision(0)) || !std::isfinite(norm)) {
                        throw std::runtime_error("DSLOAD TRVEC direction must be finite and nonzero");
                    }
                    unit_direction /= norm;

                    model.add_dload(
                        surface,
                        scale * magnitude * unit_direction,
                        *orientation,
                        stored_amplitude
                    );
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
