/**
 * @file register_dload.inl
 * @brief Registers Abaqus gravity loads.
 *
 * The supported Abaqus `DLOAD, GRAV` form combines a scalar acceleration with
 * its three-component direction and creates a FEMaster volume load on an
 * element or element set. Procedure-specific amplitude handling determines
 * whether scaling is sampled immediately or retained for dynamic evaluation.
 *
 * Other distributed-load labels and imaginary loading are outside the current
 * translation and produce explicit diagnostics.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <memory>
#include <string>

#include "../reference.h"
#include "../parser_abq.h"
#include "../../../bc/load_v.h"
#include "../../../loadcase/loadcase.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

inline void register_dload(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("DLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));
        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("AMPLITUDE").optional()
                .flag("REAL")
                .flag("IMAGINARY")
        );
        command.on_enter([&parser, amplitude](const fem::io::dsl::Keys& keys) {
            auto* loadcase = parser.active_loadcase();
            logging::error(parser.abaqus_state().step_active && loadcase != nullptr,
                "DLOAD: must appear after a supported procedure inside STEP");
            logging::error(loadcase->type_name() != "EIGENFREQ",
                "DLOAD: not supported in a FREQUENCY step");
            logging::error(!(keys.has("REAL") && keys.has("IMAGINARY")),
                "DLOAD: REAL and IMAGINARY are mutually exclusive");
            logging::error(!keys.has("IMAGINARY"),
                "DLOAD: IMAGINARY is not supported");
            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").on_empty(std::string{"EALL"})
                    .one<std::string>().name("TYPE")
                    .one<Precision>().name("MAGNITUDE")
                    .fixed<Precision, 3>().name("DIRECTION")
                )
                .bind([&parser, amplitude](const std::string& target,
                                           const std::string& type,
                                           Precision magnitude,
                                           const std::array<Precision, 3>& direction) {
                    logging::error(type == "GRAV",
                        "DLOAD: only GRAV is supported");

                    const auto [scale, resolved_amplitude] = parser.resolve_load_amplitude(*amplitude);
                    const Vec3 gravity = scale * magnitude * Vec3{direction[0], direction[1], direction[2]};
                    auto& model = parser.model();

                    model::ElementRegion::Ptr region;
                    if (model._data->elem_sets.has(target)) {
                        region = model._data->elem_sets.get(target);
                    } else {
                        region = std::make_shared<model::ElementRegion>("INTERNAL");
                        region->add(io::reader::compiled_element_id(model, target));
                    }

                    bc::Amplitude::Ptr amplitude_ptr = nullptr;
                    if (!resolved_amplitude.empty()) {
                        logging::error(model._data->amplitudes.has(resolved_amplitude),
                            "DLOAD: amplitude ", resolved_amplitude, " does not exist");
                        amplitude_ptr = model._data->amplitudes.get(resolved_amplitude);
                    }

                    auto load = std::make_shared<bc::VLoad>();
                    load->region_    = std::move(region);
                    load->values_    = gravity;
                    load->amplitude_ = std::move(amplitude_ptr);
                    model.add_load(std::move(load));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
