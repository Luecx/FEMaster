/**
 * @file register_dsload.inl
 * @brief Registers Abaqus pressure and traction surface loads.
 *
 * `DSLOAD` translates scalar `P` pressure and vector `TRVEC` traction entries on
 * compiled surface sets into FEMaster surface loads. It resolves optional
 * amplitudes and coordinate orientations, normalizes traction directions and
 * stores the resulting objects in the active step collector.
 *
 * Follower behavior is validated against the active procedure because the
 * nonlinear formulations currently support only the explicitly representable
 * combinations. Complex-valued loading is rejected.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>

#include "../parser_abq.h"
#include "../../../bc/load_d.h"
#include "../../../bc/load_p.h"
#include "../../../loadcase/loadcase.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands_abq {

inline void register_dsload(fem::io::dsl::Registry& registry, ParserAbq& parser) {
    registry.command("DSLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("STEP"));

        auto amplitude   = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();
        auto follower    = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("AMPLITUDE").optional()
                .key("ORIENTATION").optional()
                .key("FOLLOWER").optional("YES").allowed({"YES", "NO"})
                .flag("REAL")
                .flag("IMAGINARY")
        );
        command.on_enter([&parser, amplitude, orientation, follower](const fem::io::dsl::Keys& keys) {
            auto* loadcase = parser.active_loadcase();
            logging::error(parser.abaqus_state().step_active && loadcase != nullptr,
                "DSLOAD: must appear after a supported procedure inside STEP");
            logging::error(loadcase->type_name() != "EIGENFREQ",
                "DSLOAD: not supported in a FREQUENCY step");
            logging::error(!(keys.has("REAL") && keys.has("IMAGINARY")),
                "DSLOAD: REAL and IMAGINARY are mutually exclusive");
            logging::error(!keys.has("IMAGINARY"),
                "DSLOAD: IMAGINARY is not supported");

            *amplitude   = keys.has("AMPLITUDE")   ? keys.raw("AMPLITUDE")   : std::string{};
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            *follower    = keys.raw("FOLLOWER");

            logging::error(orientation->empty() || parser.model()._data->coordinate_systems.has(*orientation),
                "DSLOAD: orientation ", *orientation, " does not exist");
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("SURFACE")
                    .one<std::string>().name("TYPE")
                    .one<Precision>().name("MAGNITUDE")
                    .fixed<Precision, 3>().name("DIRECTION")
                        .on_missing(std::numeric_limits<Precision>::quiet_NaN())
                        .on_empty(std::numeric_limits<Precision>::quiet_NaN())
                )
                .bind([&parser, amplitude, orientation, follower](const std::string& surface,
                                                                  const std::string& type,
                                                                  Precision magnitude,
                                                                  const std::array<Precision, 3>& direction) {
                    auto& model = parser.model();
                    logging::error(model._data->surface_sets.has(surface),
                        "DSLOAD: surface ", surface, " does not exist");

                    const auto [scale, resolved_amplitude] = parser.resolve_load_amplitude(*amplitude);
                    magnitude *= scale;
                    const std::string procedure = parser.active_loadcase()->type_name();

                    bc::Amplitude::Ptr amplitude_ptr = nullptr;
                    if (!resolved_amplitude.empty()) {
                        logging::error(model._data->amplitudes.has(resolved_amplitude),
                            "DSLOAD: amplitude ", resolved_amplitude, " does not exist");
                        amplitude_ptr = model._data->amplitudes.get(resolved_amplitude);
                    }

                    if (type == "P") {
                        logging::error(std::isnan(direction[0]) && std::isnan(direction[1]) && std::isnan(direction[2]),
                            "DSLOAD: P accepts no direction components");
                        logging::error(procedure != "NONLINEARSTATIC",
                            "DSLOAD: follower pressure is not supported in nonlinear steps");

                        auto load = std::make_shared<bc::PLoad>();
                        load->region_    = model._data->surface_sets.get(surface);
                        load->pressure_  = magnitude;
                        load->amplitude_ = std::move(amplitude_ptr);
                        model.add_load(std::move(load));
                        return;
                    }

                    logging::error(type == "TRVEC",
                        "DSLOAD: only P and TRVEC are supported");
                    logging::error(!std::isnan(direction[0]) && !std::isnan(direction[1]) && !std::isnan(direction[2]),
                        "DSLOAD: TRVEC requires three direction components");
                    logging::error(procedure != "NONLINEARSTATIC" || *follower == "NO",
                        "DSLOAD: nonlinear TRVEC requires FOLLOWER=NO");

                    Vec3 traction{direction[0], direction[1], direction[2]};
                    const Precision norm = traction.norm();
                    logging::error(norm > Precision(0) && std::isfinite(norm),
                        "DSLOAD: TRVEC direction must be finite and nonzero");
                    traction *= magnitude / norm;

                    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
                    if (!orientation->empty()) {
                        orientation_ptr = model._data->coordinate_systems.get(*orientation);
                    }

                    auto load = std::make_shared<bc::DLoad>();
                    load->region_      = model._data->surface_sets.get(surface);
                    load->values_      = traction;
                    load->orientation_ = std::move(orientation_ptr);
                    load->amplitude_   = std::move(amplitude_ptr);
                    model.add_load(std::move(load));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands_abq
