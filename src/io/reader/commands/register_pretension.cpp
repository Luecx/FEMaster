#include "register_functions.h"
#include "../parser.h"
#include "../../dsl/registry.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"
#include "../../../bc/structural/pretension.h"

#include <memory>
#include <string>

namespace fem::io::reader::commands {

void register_pretension_section(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("PRETENSIONSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Pair two coincident, topologically separate solid cut surfaces.");
        command.keyword(fem::io::dsl::KeywordSpec::make()
            .key("NAME").required()
            .key("SURFACE_A").required()
            .key("SURFACE_B").required());
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            model.add_pretension_interface_section(
                keys.raw("NAME"),
                keys.raw("SURFACE_A"),
                keys.raw("SURFACE_B"));
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

void register_pretension(fem::io::dsl::Registry& registry, model::Model& model, Parser& parser) {
    registry.command("PRETENSION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "LOADCASE", "STATIC"}));
        command.doc("Set FORCE or HOLD on a pretension section in a linear static step.");
        command.keyword(fem::io::dsl::KeywordSpec::make()
            .key("SECTION").required()
            .key("MODE").required().allowed({"FORCE", "HOLD"})
            .key("VALUE").optional(""));
        command.on_enter([&model, &parser](const fem::io::dsl::Keys& keys) {
            const std::string name = keys.raw("SECTION");
            pretension::PretensionSection::Ptr section;
            for (const auto& candidate : model._data->pretension_sections)
                if (candidate->name == name) section = candidate;
            logging::error(section != nullptr, "PRETENSION: unknown section ", name);
            using Prescription = bc::PretensionCondition;
            Prescription::Mode mode = Prescription::Mode::Hold;
            Precision value = Precision(0);
            if (keys.raw("MODE") == "FORCE") {
                logging::error(!keys.raw("VALUE").empty(), "PRETENSION FORCE requires VALUE");
                mode = Prescription::Mode::Force;
                value = keys.get<Precision>("VALUE");
            } else {
                logging::error(keys.raw("VALUE").empty(), "PRETENSION HOLD does not accept VALUE");
            }
            parser.modify_conditions(bc::PRETENSION, name,
                {std::make_shared<Prescription>(section, mode, value)});
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
