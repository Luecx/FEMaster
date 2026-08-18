/**
 * @file register_dload.inl
 * @brief Registers supported Abaqus *DLOAD element-based distributed loads.
 *
 * The supported DLOAD subset contains gravity loading (`GRAV`). Each definition
 * is translated directly into a FEMaster volume load in the load collector of the
 * single supported Abaqus analysis step.
 *
 * @see ParserAbqState
 * @see model::Model::add_vload
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <array>
#include <charconv>
#include <memory>
#include <string>
#include <system_error>

#include "../parser_abq.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers the supported Abaqus DLOAD syntax.
 *
 * Each GRAV record contains an element or element-set target, a scalar magnitude
 * and three direction components. A blank target maps to `EALL`. Named amplitudes
 * and the step-level default amplitude are resolved according to the active
 * procedure. Imaginary harmonic load components are unsupported.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser containing the active analysis step.
 */
inline void register_dload(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("DLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Define Abaqus GRAV loads in the active analysis step.");

        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("AMPLITUDE").optional().doc("Optional named Abaqus amplitude")
                .flag("REAL").doc("Real harmonic load component")
                .flag("IMAGINARY").doc("Unsupported imaginary harmonic load component")
        );

        command.on_enter([&parser, amplitude](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            logging::error(state.step_active && parser.active_loadcase(),
                "DLOAD must appear after a supported procedure inside STEP");
            logging::error(state.procedure != "EIGENFREQ",
                "DLOAD is not supported in a FREQUENCY step");
            logging::error(!(keys.has("REAL") && keys.has("IMAGINARY")),
                "DLOAD REAL and IMAGINARY are mutually exclusive");
            logging::error(!keys.has("IMAGINARY"),
                "DLOAD IMAGINARY is not supported by the real-load harmonic solver");

            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
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
                    logging::error(type == "GRAV",
                        "DLOAD currently supports only GRAV; use SURFACE + DSLOAD for surface pressure/traction");

                    const auto [scale, resolved_amplitude] = parser.resolve_load_amplitude(*amplitude);
                    const fem::Vec3 gravity = scale * magnitude * fem::Vec3{
                        direction[0], direction[1], direction[2]
                    };

                    auto& model = parser.model();
                    if (model._data->elem_sets.has(target)) {
                        model.add_vload(target, gravity, "", resolved_amplitude);
                        return;
                    }

                    fem::ID element_id{};
                    const char* begin = target.data();
                    const char* end   = begin + target.size();
                    const auto [ptr, ec] = std::from_chars(begin, end, element_id);
                    logging::error(ec == std::errc{} && ptr == end,
                        "DLOAD target '", target, "' is not an element set or element id");
                    model.add_vload(element_id, gravity, "", resolved_amplitude);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
