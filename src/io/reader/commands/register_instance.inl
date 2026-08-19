/**
 * @file register_instance.inl
 * @brief Registers rigid part instances for FEMaster decks.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <cmath>

#include "../../../cos/rectangular_system.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

inline void register_instance(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("INSTANCE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
        command.doc("Create a rigid instance of a previously defined part.");
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").required().doc("Unique instance name")
                .key("PART").required().doc("Referenced part name")
                .key("TX").optional("0").doc("Global X translation")
                .key("TY").optional("0").doc("Global Y translation")
                .key("TZ").optional("0").doc("Global Z translation")
                .key("RX").optional("0").doc("Extrinsic X rotation in degrees")
                .key("RY").optional("0").doc("Extrinsic Y rotation in degrees")
                .key("RZ").optional("0").doc("Extrinsic Z rotation in degrees")
        );
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            const Precision deg = Precision(std::acos(-1.0) / 180.0);
            const Vec3 translation{
                keys.get<Precision>("TX"),
                keys.get<Precision>("TY"),
                keys.get<Precision>("TZ")
            };
            const Precision rx = keys.get<Precision>("RX") * deg;
            const Precision ry = keys.get<Precision>("RY") * deg;
            const Precision rz = keys.get<Precision>("RZ") * deg;
            const Mat3 rotation =
                cos::RectangularSystem::rotation_z(rz)
              * cos::RectangularSystem::rotation_y(ry)
              * cos::RectangularSystem::rotation_x(rx);

            model.add_instance(keys.raw("NAME"), keys.raw("PART"), translation, rotation);
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
