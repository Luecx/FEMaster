// Native face-pair pretension definitions and linear-static step actions.

#include "register_functions.h"
#include "../parser.h"
#include "../../dsl/registry.h"

#include <memory>
#include <string>

#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../core/types_num.h"
#include "../../../model/model.h"

namespace fem::io::reader { class Parser; }

namespace fem::io::reader::commands {

void register_pretension_section(
    fem::io::dsl::Registry& registry,
    model::Model& model) {
    registry.command("PRETENSIONSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Define a pretension section using one merged interface face "
            "or two existing face sets.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME")
                    .required()
                    .doc("Pretension section name")
                .key("SURFACE_A")
                    .required()
                    .doc("First existing interface face set")
                .key("NODE").optional().doc("Optional existing scalar control node or singleton NSET")
                .key("SURFACE_B")
                    .optional("")
                    .doc("Second interface face set; omit for a merged interface"));

        auto name = std::make_shared<std::string>();
        auto surface_a = std::make_shared<std::string>();
        auto node = std::make_shared<std::string>();
        auto surface_b = std::make_shared<std::string>();

        command.on_enter([name, surface_a, surface_b, node](const fem::io::dsl::Keys& keys) {
            *node = keys.raw("NODE");
            *name = keys.raw("NAME");
            *surface_a = keys.raw("SURFACE_A");
            *surface_b = keys.raw("SURFACE_B");
            if (surface_b->empty()) *surface_b = *surface_a;
        });

        command.on_exit([&model, name, surface_a, surface_b, node](const fem::io::dsl::Keys&) {
            model.add_pretension_interface_section(*name, *surface_a, *surface_b);
            if (!node->empty()) {
                auto region = model.resolve_node_region(*node);
                logging::error(region->size() == 1,
                    "PRETENSION SECTION: NODE must resolve to exactly one node");
                model._data->pretension_sections.back()->control_node = *region->begin();
            }
        });

        command.variant(fem::io::dsl::Variant::make());
    });
    registry.command("PRE-TENSIONSECTION", [&](fem::io::dsl::Command& command) {
        namespace dsl = fem::io::dsl;
        command.allow_if(dsl::Condition::parent_is({"ROOT", "ASSEMBLY"}));
        command.keyword(dsl::KeywordSpec::make()
            .key("SURFACE").optional()
            .key("ELEMENT").optional()
            .key("NODE").required()
            .key("FOLLOWER").optional("NO"));
        auto surface = std::make_shared<std::string>();
        auto node = std::make_shared<std::string>();
        auto normal = std::make_shared<Vec3>(Vec3::Zero());
        command.on_enter([surface, node, normal](const dsl::Keys& keys) {
            logging::error(!keys.has("ELEMENT"),
                "PRE-TENSION SECTION: ELEMENT-based pretension is unsupported; use SURFACE");
            logging::error(keys.has("SURFACE"), "PRE-TENSION SECTION: SURFACE is required");
            logging::error(keys.raw("FOLLOWER") == "NO",
                "PRE-TENSION SECTION: only FOLLOWER=NO is supported");
            *surface = keys.raw("SURFACE");
            *node = keys.raw("NODE");
            normal->setZero();
        });
        command.variant(dsl::Variant::make().segment(dsl::Segment::make()
            .range(dsl::LineRange{}.min(0).max(1))
            .pattern(dsl::Pattern::make().one<Precision>().one<Precision>().one<Precision>())
            .bind([normal](Precision x, Precision y, Precision z) {
                *normal = Vec3{x,y,z};
                logging::error(normal->allFinite() && normal->norm() > Precision(0),
                    "PRE-TENSION SECTION: normal must be finite and nonzero");
            })));
        command.on_exit([&model, surface, node, normal](const dsl::Keys&) {
            // Require an explicitly defined node: avoids modifying frozen instance maps.
            auto region = model.resolve_node_region(*node);
            logging::error(region->size() == 1,
                "PRE-TENSION SECTION: NODE must resolve to exactly one existing node");
            model.add_pretension_interface_section("__ABAQUS_PT_" + *node, *surface, *surface);
            auto& section = *model._data->pretension_sections.back();
            section.control_node = *region->begin();
            section.nodal_control = true;
            section.normal_override = *normal;
        });
    });
}

void register_pretension(
    fem::io::dsl::Registry& registry,
    model::Model& model,
    fem::io::reader::Parser& parser) {
    registry.command("PRETENSION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is({"ROOT", "LOADCASE"}));
        command.doc("Apply or lock a pretension section.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("SECTION")
                    .required()
                    .doc("Pretension section name")
                .key("ACTION")
                    .optional("LOCK")
                    .allowed({"LOCK", "LOAD"})
                .key("CONTROL")
                    .optional("DISPLACEMENT")
                    .allowed({"FORCE", "DISPLACEMENT"})
                .key("VALUE")
                    .optional("")
                    .doc("Force or relative displacement value"));

        command.on_enter([&model, &parser](const fem::io::dsl::Keys& keys) {
            const std::string section = keys.raw("SECTION");
            const std::string action = keys.raw("ACTION");
            logging::error(action != "LOAD" || !keys.raw("VALUE").empty(),
                "PRETENSION ACTION=LOAD requires VALUE");

            if (parser.active_loadcase()) {
                parser.queue_pretension_action(
                    section, action, keys.raw("CONTROL"),
                    action == "LOAD" ? keys.get<fem::Precision>("VALUE") : Precision(0));
                return;
            }

            if (action == "LOCK") {
                model.lock_pretension_section(section);
                return;
            }

            const std::string control = keys.raw("CONTROL");
            const fem::Precision value = keys.get<fem::Precision>("VALUE");
            const auto mode = control == "FORCE"
                ? fem::pretension::Control::Force
                : fem::pretension::Control::Displacement;

            model.set_pretension_load(section, mode, value);
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
