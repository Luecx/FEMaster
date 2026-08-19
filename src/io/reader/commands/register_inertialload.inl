/**
 * @file register_inertialload.inl
 * @brief Registers rigid-body inertial loads.
 *
 * `INERTIALOAD` defines translational, centrifugal and angular-acceleration
 * contributions about a supplied center for an element region. The command also
 * records whether concentrated point masses participate in the inertia
 * calculation and stores the resulting `bc::InertialLoad` by collector.
 *
 * Mass-property evaluation and conversion of rigid-body accelerations into
 * consistent element and nodal forces occur during load assembly.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <memory>
#include <string>

#include "../../../bc/load_inertial.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_inertialload(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("INERTIALOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));

        auto consider_point_masses = std::make_shared<bool>(false);
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR").required()
                .key("CONSIDER_POINT_MASSES").optional("0")
        );
        command.on_enter([&model, consider_point_masses](const fem::io::dsl::Keys& keys) {
            model._data->load_cols.activate(keys.raw("LOAD_COLLECTOR"));
            *consider_point_masses = keys.get<bool>("CONSIDER_POINT_MASSES");
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .fixed<fem::Precision, 3>().name("CENTER")
                    .fixed<fem::Precision, 3>().name("CENTER_ACC")
                    .fixed<fem::Precision, 3>().name("OMEGA")
                    .fixed<fem::Precision, 3>().name("ALPHA")
                )
                .bind([&model, consider_point_masses](const std::string& target,
                                                       const std::array<fem::Precision, 3>& center,
                                                       const std::array<fem::Precision, 3>& center_acc,
                                                       const std::array<fem::Precision, 3>& omega,
                                                       const std::array<fem::Precision, 3>& alpha) {
                    logging::error(model._data->elem_sets.has(target),
                        "INERTIALOAD: element set ", target, " does not exist");

                    auto load = std::make_shared<bc::InertialLoad>();
                    load->region_                = model._data->elem_sets.get(target);
                    load->center_                = Vec3{center[0], center[1], center[2]};
                    load->center_acc_            = Vec3{center_acc[0], center_acc[1], center_acc[2]};
                    load->omega_                 = Vec3{omega[0], omega[1], omega[2]};
                    load->alpha_                 = Vec3{alpha[0], alpha[1], alpha[2]};
                    load->consider_point_masses_ = *consider_point_masses;
                    model.add_load(std::move(load));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
