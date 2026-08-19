/**
 * @file register_material.inl
 * @brief Registers named material objects for FEMaster input decks.
 *
 * The parser constructs the concrete `material::Material` container with its
 * intrinsic name and registers the finished object through
 * `Model::add_material()`. Subsequent child material-property commands operate
 * on the dictionary's current material and attach elasticity, density, thermal
 * expansion or other constitutive components directly to that object.
 *
 * @see material::Material
 * @see model::Model::add_material
 *
 * @author Finn Eggers
 * @date 18.08.2026
 */

#include <memory>
#include <string>

#include "../../../material/material.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

inline void register_material(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("MATERIAL", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Create and activate a named material definition.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("MATERIAL")
                    .alternative("NAME")
                    .required()
                    .doc("Identifier of the material entry")
        );

        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            model.add_material(std::make_shared<material::Material>(keys.raw("MATERIAL")));
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands