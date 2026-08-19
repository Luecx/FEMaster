/**
 * @file register_instance.inl
 * @brief Registers Abaqus part-instance placement syntax.
 *
 * Abaqus `INSTANCE` references a semantic Part and records a rigid translation
 * plus optional axis-angle rotation. The transformation is represented by a
 * FEMaster rectangular coordinate system attached to the new Instance rather
 * than applied destructively to part-local coordinates.
 *
 * The compile step later expands the reusable Part, transforms its nodes and
 * builds maps used by qualified assembly references.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <memory>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

inline void register_instance(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("INSTANCE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
        command.doc("Create a rigid Abaqus instance of a previously defined part.");

        auto instance = std::make_shared<model::Instance::Ptr>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").required().doc("Unique instance name")
                .key("PART").required().doc("Referenced part name")
        );

        command.on_enter([&model, instance](const fem::io::dsl::Keys& keys) {
            model.add_instance(keys.raw("NAME"), keys.raw("PART"));
            *instance = model._data->instances.get();
            logging::error(*instance != nullptr,
                "INSTANCE: active instance was not created");
        });

        auto rotation = [](const Vec3& point_a, const Vec3& point_b, Precision angle) {
            const Vec3 axis = point_b - point_a;
            logging::error(axis.norm() > std::numeric_limits<Precision>::epsilon(),
                "INSTANCE: rotation axis points must be distinct");

            const Vec3 unit = axis.normalized();
            Mat3 skew;
            skew << Precision(0), -unit(2), unit(1),
                    unit(2), Precision(0), -unit(0),
                   -unit(1), unit(0), Precision(0);

            const Precision c = std::cos(angle);
            const Precision s = std::sin(angle);
            return c * Mat3::Identity()
                 + (Precision(1) - c) * (unit * unit.transpose())
                 + s * skew;
        };

        command.variant(fem::io::dsl::Variant::make()
            .rank(30)
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(2))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 10>().name("PLACEMENT")
                )
                .bind([instance, rotation](const std::array<fem::Precision, 10>& values) {
                    logging::error(*instance != nullptr,
                        "INSTANCE: internal instance state is not initialized");

                    const Vec3 translation{values[0], values[1], values[2]};
                    const Vec3 point_a{values[3], values[4], values[5]};
                    const Vec3 point_b{values[6], values[7], values[8]};
                    const Precision angle = values[9] * Precision(std::acos(-1.0) / 180.0);
                    const Mat3 R = rotation(point_a, point_b, angle);

                    (*instance)->rotation    = R;
                    (*instance)->translation = R * translation + point_a - R * point_a;
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .rank(20)
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 7>().name("ROTATION")
                )
                .bind([instance, rotation](const std::array<fem::Precision, 7>& values) {
                    logging::error(*instance != nullptr,
                        "INSTANCE: internal instance state is not initialized");

                    const Vec3 point_a{values[0], values[1], values[2]};
                    const Vec3 point_b{values[3], values[4], values[5]};
                    const Precision angle = values[6] * Precision(std::acos(-1.0) / 180.0);
                    const Mat3 R = rotation(point_a, point_b, angle);

                    (*instance)->rotation    = R;
                    (*instance)->translation = point_a - R * point_a;
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .rank(10)
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 3>().name("TRANSLATION")
                )
                .bind([instance](const std::array<fem::Precision, 3>& values) {
                    logging::error(*instance != nullptr,
                        "INSTANCE: internal instance state is not initialized");
                    (*instance)->translation = Vec3{values[0], values[1], values[2]};
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make().rank(0));
    });
}

} // namespace fem::io::reader::commands_abq
