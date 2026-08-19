/**
 * @file register_support.inl
 * @brief Registers nodal supports.
 *
 * A `SUPPORT` row targets a node or node set and supplies six generalized
 * constraint values. Finite entries prescribe translations or rotations while
 * omitted components remain unconstrained; an optional coordinate system
 * defines the basis in which the values are expressed.
 *
 * The completed `bc::Support` is stored in the named support collector. Actual
 * constraint equations are assembled later by the active load case.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <limits>
#include <memory>
#include <string>

#include "../reference.h"
#include "../../../bc/support.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

inline void register_support(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("SUPPORT", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));

        auto orientation = std::make_shared<std::string>();
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("SUPPORT_COLLECTOR").required()
                .key("ORIENTATION").optional()
        );
        command.on_enter([&model, orientation](const fem::io::dsl::Keys& keys) {
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            model._data->supp_cols.activate(keys.raw("SUPPORT_COLLECTOR"));
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET")
                    .fixed<fem::Precision, 6>().name("DOF")
                        .on_missing(std::numeric_limits<fem::Precision>::quiet_NaN())
                        .on_empty(std::numeric_limits<fem::Precision>::quiet_NaN())
                )
                .bind([&model, orientation](const std::string& target,
                                            const std::array<fem::Precision, 6>& values) {
                    model::NodeRegion::Ptr region;
                    if (model._data->node_sets.has(target)) {
                        region = model._data->node_sets.get(target);
                    } else {
                        region = std::make_shared<model::NodeRegion>("INTERNAL");
                        region->add(io::reader::compiled_node_id(model, target));
                    }

                    cos::CoordinateSystem::Ptr orientation_ptr = nullptr;
                    if (!orientation->empty()) {
                        logging::error(model._data->coordinate_systems.has(*orientation),
                            "SUPPORT: coordinate system ", *orientation, " does not exist");
                        orientation_ptr = model._data->coordinate_systems.get(*orientation);
                    }

                    Vec6 constraint;
                    for (Index i = 0; i < 6; ++i) {
                        constraint(i) = values[static_cast<std::size_t>(i)];
                    }
                    model.add_support(bc::Support{std::move(region), constraint, std::move(orientation_ptr)});
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
