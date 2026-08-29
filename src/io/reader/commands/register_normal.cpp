/**
 * @file register_normal.cpp
 * @brief Registers shell-normal field selection.
 *
 * The `*NORMAL` command selects an existing three-component
 * `ELEMENT_NODAL` field and publishes it through
 * `ModelData::shell_element_nodal_normals`.
 *
 * The command does not create, populate, normalize or average normal values.
 * Geometric completion and angular equalisation remain responsibilities of
 * `Model::build_shell_element_normals()`.
 *
 * @see model::Field
 * @see model::Model::build_shell_element_normals
 *
 * @author Finn Eggers
 * @date 30.07.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../../../core/logging.h"
#include "../../../data/field.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

#include <string>

namespace fem::io::reader::commands {

/**
 * Registers selection of the element-nodal field used for shell normals.
 *
 * `*NORMAL, FIELD=<name>` resolves one previously created generic field,
 * validates the required domain and component count and directly assigns the
 * corresponding ModelData pointer. Field completion is performed later by the
 * existing model-side shell-normal build routine.
 *
 * @param registry DSL registry receiving the `*NORMAL` command definition.
 * @param model Model whose shell-normal field pointer is modified.
 */
void register_normal(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("NORMAL", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Use an ELEMENT_NODAL vector field as the shell normal field.");
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("FIELD").required().doc("ELEMENT_NODAL field containing shell normals")
        );

        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            // Resolve the already populated generic field. The command only
            // selects this storage; completion remains model-side.
            const std::string field_name = keys.raw("FIELD");
            const auto        field      = model._data->get_field(field_name);

            logging::error(field != nullptr,
                "NORMAL: field '", field_name, "' does not exist");
            logging::error(field->domain == model::FieldDomain::ELEMENT_NODAL,
                "NORMAL: field '", field_name, "' must use ELEMENT_NODAL domain");
            logging::error(field->components == 3,
                "NORMAL: field '", field_name, "' must have exactly three components");

            // Publish the field directly. Existing entries are interpreted by
            // Model::build_shell_element_normals() after the field parser stage.
            model._data->shell_element_nodal_normals = field;
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
