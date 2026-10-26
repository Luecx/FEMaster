#include "register_functions.h"
#include "../parser.h"
#include "../../dsl/registry.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

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
        auto name = std::make_shared<std::string>();
        auto surface_a = std::make_shared<std::string>();
        auto surface_b = std::make_shared<std::string>();
        command.on_enter([name, surface_a, surface_b](const fem::io::dsl::Keys& keys) {
            *name = keys.raw("NAME");
            *surface_a = keys.raw("SURFACE_A");
            *surface_b = keys.raw("SURFACE_B");
        });
        command.on_exit([&model, name, surface_a, surface_b](const fem::io::dsl::Keys&) {
            model.add_pretension_interface_section(*name, *surface_a, *surface_b);
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

void register_pretension(fem::io::dsl::Registry& registry, model::Model& model, Parser&) {
    registry.command("PRETENSION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "LOADCASE", "STATIC"}));
        command.doc("Set FORCE or HOLD on a pretension section in a linear static step.");
        command.keyword(fem::io::dsl::KeywordSpec::make()
            .key("SECTION").required()
            .key("MODE").required().allowed({"FORCE", "HOLD"})
            .key("VALUE").optional(""));
        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            const std::string section = keys.raw("SECTION");
            if (keys.raw("MODE") == "FORCE") {
                logging::error(!keys.raw("VALUE").empty(), "PRETENSION FORCE requires VALUE");
                model.set_pretension_load(section, keys.get<Precision>("VALUE"));
            } else {
                logging::error(keys.raw("VALUE").empty(), "PRETENSION HOLD does not accept VALUE");
                model.lock_pretension_section(section);
            }
        });
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands