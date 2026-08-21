/**
 * @file register_instance.inl
 * @brief Registers rigid part instances for FEMaster and Abaqus decks.
 *
 * `INSTANCE` references a previously defined Part. Its optional first data line
 * contains the translation vector, and its optional second data line describes
 * a rotation axis through two global points followed by an angle in degrees.
 *
 * The parser composes both operations into the affine rigid placement stored by
 * `model::Instance`:
 *
 * \f[
 *     \boldsymbol x = \boldsymbol R\,\boldsymbol X
 *                   + \boldsymbol R\,\boldsymbol t
 *                   + \boldsymbol p-\boldsymbol R\,\boldsymbol p.
 * \f]
 *
 * The source Part remains unchanged. `Model::compile()` later expands every
 * instance, applies the stored rotation and effective translation and constructs
 * deterministic maps from local to dense assembly identifiers.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <memory>

#include "../../../core/logging.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands_abq {

/**
 * Registers the shared FEMaster/Abaqus instance-placement grammar.
 *
 * The command creates an identity placement when it is entered and then selects
 * one of four data layouts: no placement, translation only, rotation only or a
 * translation line followed by a rotation line. Rotation axes are specified by
 * two points in global assembly coordinates and converted to a proper rotation
 * with Rodrigues' formula.
 *
 * For the combined form, the translation acts first. Rotation about an axis
 * through `point_a` then produces the effective affine translation
 * `R * translation + point_a - R * point_a`, which is stored together with `R`
 * on the semantic Instance. No part-local coordinate is modified while parsing.
 *
 * @param registry Parser registry receiving the `INSTANCE` command definition.
 * @param model Model owning Parts and the newly constructed Instance.
 */
inline void register_instance(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("INSTANCE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
        command.doc(
            "Create a rigid instance of a previously defined part. The optional "
            "first data line contains tx, ty, tz. The optional second line contains "
            "x1, y1, z1, x2, y2, z2, angle in degrees."
        );

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

        auto build_rotation = [](const Vec3& point_a, const Vec3& point_b, Precision angle) {
            // Validate and normalize the global rotation axis
            logging::error(point_a.allFinite() && point_b.allFinite() && std::isfinite(angle),
                "INSTANCE: rotation definition must contain only finite values");

            const Vec3 axis = point_b - point_a;
            logging::error(axis.norm() > std::numeric_limits<Precision>::epsilon(),
                "INSTANCE: rotation axis points must be distinct");
            const Vec3 unit = axis.normalized();

            // Evaluate every component of Rodrigues' formula explicitly. This
            // keeps the fixed-size Eigen result fully initialized without an
            // intermediate expression involving a skew-symmetric matrix.
            const Precision x = unit(0);
            const Precision y = unit(1);
            const Precision z = unit(2);
            const Precision c = std::cos(angle);
            const Precision s = std::sin(angle);
            const Precision t = Precision(1) - c;

            Mat3 rotation;
            rotation(0, 0) = c + t * x * x;
            rotation(0, 1) = t * x * y - s * z;
            rotation(0, 2) = t * x * z + s * y;
            rotation(1, 0) = t * x * y + s * z;
            rotation(1, 1) = c + t * y * y;
            rotation(1, 2) = t * y * z - s * x;
            rotation(2, 0) = t * x * z - s * y;
            rotation(2, 1) = t * y * z + s * x;
            rotation(2, 2) = c + t * z * z;
            return rotation;
        };

        // Translation followed by rotation about a global two-point axis
        command.variant(fem::io::dsl::Variant::make()
            .rank(30)
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
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 7>().name("ROTATION")
                )
                .bind([instance, build_rotation](const std::array<fem::Precision, 7>& values) {
                    logging::error(*instance != nullptr,
                        "INSTANCE: internal instance state is not initialized");

                    // Convert the Abaqus axis-angle definition to affine R and t
                    const Vec3 point_a{values[0], values[1], values[2]};
                    const Vec3 point_b{values[3], values[4], values[5]};
                    const Precision angle = values[6] * Precision(std::acos(-1.0) / 180.0);
                    const Mat3 rotation = build_rotation(point_a, point_b, angle);

                    (*instance)->rotation     = rotation;
                    (*instance)->translation = rotation * (*instance)->translation
                                             + point_a - rotation * point_a;
                })
            )
        );

        // Rotation without a preceding translation
        command.variant(fem::io::dsl::Variant::make()
            .rank(20)
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .fixed<fem::Precision, 7>().name("ROTATION")
                )
                .bind([instance, build_rotation](const std::array<fem::Precision, 7>& values) {
                    logging::error(*instance != nullptr,
                        "INSTANCE: internal instance state is not initialized");

                    // Include the off-origin axis point in the affine translation
                    const Vec3 point_a{values[0], values[1], values[2]};
                    const Vec3 point_b{values[3], values[4], values[5]};
                    const Precision angle = values[6] * Precision(std::acos(-1.0) / 180.0);
                    const Mat3 rotation = build_rotation(point_a, point_b, angle);

                    (*instance)->rotation    = rotation;
                    (*instance)->translation = point_a - rotation * point_a;
                })
            )
        );

        // Translation without a following rotation
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

        // Identity placement when the instance has no data lines
        command.variant(fem::io::dsl::Variant::make().rank(0));
    });
}

} // namespace fem::io::reader::commands_abq
