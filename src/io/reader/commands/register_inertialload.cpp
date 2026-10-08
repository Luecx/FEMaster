/**
 * @file register_inertialload.cpp
 * @brief Registers rigid-body inertial loads.
 *
 * `INERTIALOAD` defines translational, centrifugal and angular-acceleration
 * contributions about a supplied center for a compiled element region. The command
 * also records whether concentrated point masses participate in the inertia
 * calculation. Explicit NAME definitions enter reusable collector storage;
 * unnamed analysis rows apply only during their active step.
 * Model resolves named element sets and scalar references into compiled regions.
 *
 * Mass-property evaluation and conversion of rigid-body accelerations into
 * consistent element and nodal forces occur during load assembly.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <array>
#include <memory>
#include <string>
#include <utility>

#include "../parser.h"
#include "../../../bc/structural/load_inertial.h"
#include "../../../core/logging.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * @brief Registers rigid-body acceleration fields on compiled element regions.
 *
 * NAME selects an explicit collector outside LOADCASE or an optional collector
 * inside an active analysis. Without NAME, the load is step-local, while
 * OP=NEW retires inherited INERTIAL_LOAD history before reading the block.
 * Targets resolve to named sets or scalar references.
 * Center, translational acceleration, angular velocity and angular acceleration
 * retain their global-coordinate convention; nodal TRANSFORM is not applied.
 * The created load evaluates mass and consistent inertial forces during assembly.
 *
 * @param registry Registry receiving the grammar and callbacks.
 * @param parser Parser supplying the compiled model and active analysis.
 */
void register_inertialload(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("INERTIALOAD", [&](fem::io::dsl::Command& command) {
        // Inertial loads require compiled elements and may belong to an active analysis.
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY", "LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        command.doc(
            "Create FEMaster rigid-body inertial loads; Abaqus-compatible gravity "
            "is available through DLOAD, GRAV."
        );

        auto consider_point_masses = std::make_shared<bool>(false);
        // An explicit NAME selects reusable definition storage. Unnamed analysis
        // rows update direct history and never inherit a prior insertion target.
        auto collector = std::make_shared<bc::LoadCollector::Ptr>(nullptr);

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("OP").optional("MOD").allowed({"MOD", "NEW"})
                .key("NAME").alternative("LOAD_COLLECTOR").alternative("LOADCOLLECTOR")
                    .optional().doc("Collector name; optional inside LOADCASE")
                .key("CONSIDER_POINT_MASSES").optional("0")
        );
        command.on_enter([&parser, consider_point_masses, collector](const fem::io::dsl::Keys& keys) {
            // Resolve mass participation and the definition destination once
            auto& model = parser.model();
            *consider_point_masses = keys.get<bool>("CONSIDER_POINT_MASSES");
            // NAME denotes a named definition, while an unnamed row in an
            // analysis changes only the current model history.
            const std::string collector_name = keys.raw("NAME");
            logging::error(!collector_name.empty() || parser.active_loadcase() != nullptr,
                "INERTIALLOAD: definitions outside an analysis require NAME");
            collector->reset();
            if (!collector_name.empty()) {
                parser.activate_load_collector(collector_name);
                *collector = model._data->load_cols.get();
            } else if (keys.raw("OP") == "NEW") {
                parser.clear_conditions(bc::INERTIAL_LOAD);
            }
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string        >().name("TARGET"    )
                    .fixed<fem::Precision, 3>().name("CENTER"    )
                    .fixed<fem::Precision, 3>().name("CENTER_ACC")
                    .fixed<fem::Precision, 3>().name("OMEGA"     )
                    .fixed<fem::Precision, 3>().name("ALPHA"     )
                )
                .bind([&parser, consider_point_masses, collector](
                    const std::string&                  target,
                    const std::array<fem::Precision, 3>& center,
                    const std::array<fem::Precision, 3>& center_acc,
                    const std::array<fem::Precision, 3>& omega,
                    const std::array<fem::Precision, 3>& alpha
                ) {
                    auto& model = parser.model();

                    // Keep all rigid-body kinematics in the global basis. The model
                    // maps the region while InertialLoad assembles its mass response.
                    auto load = std::make_shared<bc::InertialLoad>();
                    load->region_                = model.resolve_element_region(target);
                    load->center_                = Vec3{center    [0], center    [1], center    [2]};
                    load->center_acc_            = Vec3{center_acc[0], center_acc[1], center_acc[2]};
                    load->omega_                 = Vec3{omega     [0], omega     [1], omega     [2]};
                    load->alpha_                 = Vec3{alpha     [0], alpha     [1], alpha     [2]};
                    load->consider_point_masses_ = *consider_point_masses;

                    // Named definitions and direct history are mutually exclusive targets
                    if (*collector) {
                        (*collector)->add(std::move(load));
                    } else {
                        parser.select_collector_condition(bc::INERTIAL_LOAD, std::move(load));
                    }
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands