/**
 * @file register_dsload.cpp
 * @brief Registers Abaqus-compatible surface-based distributed loads.
 *
 * `DSLOAD` translates scalar `P` pressure entries on
 * compiled surface regions into FEMaster surface loads. It resolves optional
 * amplitudes and stores named definitions only in their explicit collector.
 * TRVEC is rejected until its semantics are implemented. Unnamed analysis
 * definitions update direct DSLOAD history with OP=MOD/NEW. Model resolves the
 * compiled surface region; the DSL owns history replacement and named selection.
 * Amplitudes follow native FEMaster evaluation during load assembly rather than
 * the Abaqus reader's procedure-dependent step scaling.
 *
 * Follower pressure is rejected in nonlinear procedures because a consistent
 * load stiffness is unavailable. Complex-valued loading is rejected.
 *
 * @author Finn Eggers
 * @date 07.10.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <utility>

#include "../parser.h"
#include "../../../bc/structural/load_d.h"
#include "../../../bc/structural/load_p.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * @brief Registers the supported Abaqus DSLOAD pressure syntax.
 *
 * P rows specify a surface and scalar pressure with no direction components.
 * TRVEC entries are explicitly rejected. ORIENTATION does not affect pressure.
 * Nodal TRANSFORM assignments do not affect these distributed surface loads.
 *
 * Entry resolves modifiers and the named definition target or direct history.
 * OP=NEW clears only direct DSLOAD history; MOD replaces matching surface targets.
 * Native amplitude semantics remain in the load objects. As in the existing
 * Abaqus registration, imaginary loads are rejected, nonlinear pressure is not
 * supported and TRVEC is rejected. This registration does not implement
 * the full Abaqus load-label or load-history parameter catalogue.
 *
 * @param registry Registry receiving the grammar and callbacks.
 * @param parser Parser supplying the compiled model and active analysis.
 */
void register_dsload(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("DSLOAD", [&](fem::io::dsl::Command& command) {
        // Match native load scope and collector semantics while retaining DSLOAD rows.
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        command.doc("Create Abaqus-compatible surface-based P pressure; TRVEC is currently unsupported.");

        auto amplitude   = std::make_shared<bc::Amplitude::Ptr>(nullptr);
        auto orientation = std::make_shared<cos::CoordinateSystem::Ptr>(nullptr);

        // An explicit NAME selects reusable definition storage. Unnamed analysis
        // rows update direct history and never inherit a prior insertion target.
        auto collector = std::make_shared<bc::LoadCollector::Ptr>(nullptr);

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
                .key("NAME").alternative("LOAD_COLLECTOR").alternative("LOADCOLLECTOR")
                    .optional().doc("Collector name; optional inside LOADCASE")
                .key("AMPLITUDE").optional()
                .key("ORIENTATION").optional()
                .key("FOLLOWER").optional("NO").allowed({"NO"})
                .flag("REAL")
                .flag("IMAGINARY")
        );
        command.on_enter([&parser, amplitude, orientation, collector](const fem::io::dsl::Keys& keys) {
            auto& model = parser.model();

            orientation->reset();
            amplitude  ->reset();
            collector  ->reset();

            const std::string orientation_name = keys.raw("ORIENTATION");
            const std::string amplitude_name   = keys.raw("AMPLITUDE");
            const std::string collector_name   = keys.raw("NAME");

            const auto* loadcase = parser.active_loadcase();
            logging::error(loadcase == nullptr || loadcase->type_name() != "EIGENFREQ",
                "DSLOAD: not supported in a FREQUENCY step");
            logging::error(!(keys.has("REAL") && keys.has("IMAGINARY")),
                "DSLOAD: REAL and IMAGINARY are mutually exclusive");
            logging::error(!keys.has("IMAGINARY"),
                "DSLOAD: IMAGINARY is not supported yet");
            logging::error(amplitude_name.empty() || model._data->amplitudes.has(amplitude_name),
                "DSLOAD: amplitude ", amplitude_name, " does not exist");
            logging::error(orientation_name.empty() || model._data->coordinate_systems.has(orientation_name),
                "DSLOAD: orientation ", orientation_name, " does not exist");
            logging::error(!collector_name.empty() || loadcase != nullptr,
                "DSLOAD: definitions outside an analysis require NAME");

            if (!amplitude_name.empty())
                *amplitude = model._data->amplitudes.get(amplitude_name);
            if (!orientation_name.empty())
                *orientation = model._data->coordinate_systems.get(orientation_name);
            if (!collector_name.empty()) {
                parser.activate_load_collector(collector_name);
                *collector = model._data->load_cols.get();
            }
            if (collector_name.empty() && keys.raw("OP") == "NEW")
                parser.clear_conditions(bc::DSLOAD);
        });

        // NaN marks omitted direction components so pressure can reject them.
        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string   >().name("SURFACE"  )
                    .one<std::string   >().name("TYPE"     )
                    .one<Precision     >().name("MAGNITUDE")
                    .fixed<Precision, 3>().name("DIRECTION")
                        .on_missing(std::numeric_limits<Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<Precision>::quiet_NaN())
                )
                .bind([&parser, amplitude, collector](
                    const std::string&              surface,
                    const std::string&              type,
                    Precision                       magnitude,
                    const std::array<Precision, 3>& direction
                ) {
                    logging::error(type != "TRVEC",
                        "DSLOAD: TRVEC is not supported");
                    logging::error(type == "P",
                        "DSLOAD: only P is supported");

                    const auto* loadcase = parser.active_loadcase();
                    logging::error(loadcase == nullptr || loadcase->type_name() != "NONLINEARSTATIC",
                        "DSLOAD: follower pressure is not supported in nonlinear steps");
                    logging::error(std::isnan(direction[0]) && std::isnan(direction[1]) && std::isnan(direction[2]),
                        "DSLOAD: P accepts no direction components");

                    // Pressure follows the surface normal. Its current formulation
                    // supports linear analyses only, including named definitions.
                    auto load = std::make_shared<bc::PLoad>();
                    load->region_    = parser.model().resolve_surface_region(surface);
                    load->pressure_  = magnitude;
                    load->amplitude_ = *amplitude;
                    if (*collector) {
                        (*collector)->add(std::move(load));
                    } else {
                        parser.modify_conditions(bc::DSLOAD, "SURFACE:" + surface + ":" + type, {std::move(load)});
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands