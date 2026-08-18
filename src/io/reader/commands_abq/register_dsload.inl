/**
 * @file register_dsload.inl
 * @brief Registers supported Abaqus *DSLOAD surface-load definitions.
 *
 * Uniform pressure (`P`) and general traction (`TRVEC`) are translated directly
 * into FEMaster surface loads in the collector of the single supported Abaqus
 * analysis step.
 *
 * @see ParserAbqState
 * @see model::Model::add_pload
 * @see model::Model::add_dload
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
 * Named amplitudes and the step-level default amplitude are resolved according to
 * the active procedure. Imaginary harmonic load components are unsupported.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser containing the active analysis step.
 */
inline void register_dsload(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("DSLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Define pressure or general traction loads in the active Abaqus step.");

        auto amplitude   = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();
        auto follower    = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("AMPLITUDE").optional().doc("Optional named Abaqus amplitude")
                .key("ORIENTATION").optional().doc("Optional traction coordinate system")
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

            *amplitude   = keys.has("AMPLITUDE")   ? keys.raw("AMPLITUDE")   : std::string{};
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            *follower    = keys.raw("FOLLOWER");

            logging::error(orientation->empty() || parser.model()._data->coordinate_systems.has(*orientation),
                "DSLOAD references unknown orientation '", *orientation, "'");
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
                .bind([&parser, amplitude, orientation, follower](const std::string& surface,
                                                                  const std::string& type,
                                                                  fem::Precision magnitude,
                                                                  const std::array<fem::Precision, 3>& direction) {
                    auto& model = parser.model();
                    auto& state = parser.abaqus_state();
                    logging::error(model._data->surface_sets.has(surface),
                        "DSLOAD references unknown surface '", surface, "'");

                    const auto [scale, resolved_amplitude] = parser.resolve_load_amplitude(*amplitude);
                    magnitude *= scale;

                    if (type == "P") {
                        logging::error(std::isnan(direction[0])
                                    && std::isnan(direction[1])
                                    && std::isnan(direction[2]),
                            "DSLOAD P accepts no traction direction components");
                        logging::error(state.procedure != "NONLINEARSTATIC"
                                    && state.procedure != "STATIC_RIKS",
                            "DSLOAD P is a follower pressure in Abaqus and is not supported in nonlinear FEMaster steps");
                        model.add_pload(surface, magnitude, resolved_amplitude);
                        return;
                    }

                    logging::error(type == "TRVEC",
                        "DSLOAD supports only P and TRVEC load labels");
                    logging::error(!std::isnan(direction[0])
                                && !std::isnan(direction[1])
                                && !std::isnan(direction[2]),
                        "DSLOAD TRVEC requires three direction components");
                    logging::error(state.procedure != "NONLINEARSTATIC"
                                && state.procedure != "STATIC_RIKS"
                                || *follower == "NO",
                        "Nonlinear DSLOAD TRVEC requires FOLLOWER=NO; follower traction is not implemented");

                    fem::Vec3 traction_direction{direction[0], direction[1], direction[2]};
                    const fem::Precision norm = traction_direction.norm();
                    logging::error(norm > fem::Precision(0) && std::isfinite(norm),
                        "DSLOAD TRVEC direction must be finite and nonzero");
                    traction_direction /= norm;

                    model.add_dload(
                        surface,
                        magnitude * traction_direction,
                        *orientation,
                        resolved_amplitude
                    );
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
