/**
 * @file register_dload.inl
 * @brief Registers the supported Abaqus *DLOAD gravity form.
 *
 * The initial element-based distributed-load reader supports `GRAV` only. The
 * Abaqus gravity magnitude and vector components map directly to FEMaster's
 * density-scaled `VLoad`, whose element integration multiplies the supplied
 * acceleration field by material density and volume measure.
 *
 * Element-face pressure/traction labels are intentionally left to `*SURFACE`
 * plus `*DSLOAD` for now. Creating new geometric surfaces during the final data
 * pass would occur after FEMaster's topology finalization and result-model write,
 * so doing that implicitly inside `*DLOAD` would violate the staged parser
 * architecture.
 *
 * @see bc::VLoad
 * @see commands_abq::register_dsload
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#pragma once

#include <array>
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
 * Registers Abaqus gravity loading on an element set, element id or the complete
 * model. The blank Abaqus target is mapped to FEMaster's built-in `EALL` region.
 *
 * @param registry Stage-local DSL registry.
 * @param parser Abaqus parser providing current step state.
 */
inline void register_dload(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("DLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        command.doc("Apply Abaqus GRAV body loading in the current step.");

        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("AMPLITUDE").optional().doc("Optional named Abaqus amplitude")
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
                    .doc("Accepted for the independent current-step collector")
        );

        command.on_enter([&parser, amplitude](const fem::io::dsl::Keys& keys) {
            auto& state = parser.abaqus_state();
            if (!state.step_active || !parser.active_loadcase()) {
                throw std::runtime_error("DLOAD must appear after a supported procedure inside STEP");
            }

            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            parser.model()._data->load_cols.activate(state.load_collector);
            state.load_collector_used = true;
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

                    auto& model = parser.model();
                    auto& state = parser.abaqus_state();
                    std::string stored_amplitude;
                    fem::Precision scale = fem::Precision(1);

                    if (!amplitude->empty()) {
                        if (!model._data->amplitudes.has(*amplitude)) {
                            throw std::runtime_error("DLOAD references unknown amplitude '" + *amplitude + "'");
                        }

                        if (state.procedure == "LINEARTRANSIENT") {
                            stored_amplitude = *amplitude;
                        } else if (state.procedure == "LINEARSTATIC" ||
                                   state.procedure == "LINEARBUCKLING") {
                            scale = model._data->amplitudes.get(*amplitude)->evaluate(state.step_period);
                        } else if (state.procedure == "NONLINEARSTATIC" ||
                                   state.procedure == "STATIC_RIKS") {
                            throw std::runtime_error(
                                "DLOAD AMPLITUDE is not supported for nonlinear static/Riks proportional loading"
                            );
                        } else if (state.procedure == "LINEARHARMONIC") {
                            throw std::runtime_error("DLOAD AMPLITUDE is not supported for harmonic response");
                        }
                    } else if (state.procedure == "LINEARTRANSIENT" &&
                               !state.step_amplitude.empty() &&
                               state.step_amplitude != "STEP") {
                        stored_amplitude = "__ABQ_STEP_" + std::to_string(state.step_index) + "_DEFAULT_AMPLITUDE";
                    }

                    const fem::Vec3 gravity = scale * magnitude * fem::Vec3{
                        direction[0], direction[1], direction[2]
                    };

                    if (model._data->elem_sets.has(target)) {
                        model.add_vload(target, gravity, "", stored_amplitude);
                        return;
                    }

                    std::size_t parsed = 0;
                    const long value = std::stol(target, &parsed);
                    if (parsed != target.size()) {
                        throw std::runtime_error("DLOAD target '" + target + "' is not an element set or element id");
                    }
                    model.add_vload(static_cast<fem::ID>(value), gravity, "", stored_amplitude);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
