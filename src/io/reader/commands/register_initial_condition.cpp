/**
 * @file register_initial_condition.cpp
 * @brief Registers model-level initial temperature and velocity fields.
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

void register_initial_condition(fem::io::dsl::Registry& registry, model::Model& model) {
    const auto register_name = [&](const char* name) {
        registry.command(name, [&](fem::io::dsl::Command& command) {
            command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
            command.doc("Assign a named FIELD as an initial model state.");

            command.keyword(
                fem::io::dsl::KeywordSpec::make()
                    .key("TYPE").required().allowed({"TEMPERATURE", "VELOCITY"})
                    .key("FIELD").required()
            );

            command.on_enter([&model](const fem::io::dsl::Keys& keys) {
                const std::string type       = keys.raw("TYPE");
                const std::string field_name = keys.raw("FIELD");
                const auto field = model._data->get_field(field_name);

                logging::error(field != nullptr,
                    "INITIALCONDITION: field ", field_name, " does not exist");
                logging::error(field->domain == model::FieldDomain::NODE,
                    "INITIALCONDITION: field ", field_name, " must use NODE domain");

                if (type == "TEMPERATURE") {
                    logging::error(field->components == 1,
                        "INITIALCONDITION: temperature field must have one component");
                    model._data->temperature = field;
                    return;
                }

                if (type == "VELOCITY") {
                    logging::error(field->components == 6,
                        "INITIALCONDITION: velocity field must have six components");
                    model._data->velocity = field;
                    return;
                }

            });

            command.variant(fem::io::dsl::Variant::make());
        });
    };

    register_name("INITIALCONDITION");
    register_name("INITIALCONDITIONS");
}

} // namespace fem::io::reader::commands
