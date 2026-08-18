// register_loadcase_topodensity.inl — registers TOPODENSITY for LINEARSTATICTOPO loadcases

#include <string>

#include "../parser.h"

#include "../../../core/logging.h"
#include "../../../loadcase/linear_static_topo.h"

namespace fem::io::reader::commands {

inline void register_loadcase_topodensity(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("TOPODENSITY", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Select the element density field for LINEARSTATICTOPO loadcases.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("FIELD").required().doc("Density field name (ELEMENT, 1 component)")
                .alternative("NAME")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto* lc = parser.active_loadcase_as<loadcase::LinearStaticTopo>();
            logging::error(lc != nullptr,
                "TOPODENSITY only valid for LINEARSTATICTOPO loadcases");

            const std::string field_name = keys.raw("FIELD");
            auto field = parser.model()._data->get_field(field_name);
            logging::error(field != nullptr,
                "TOPODENSITY field '", field_name, "' does not exist");
            logging::error(field->domain == model::FieldDomain::ELEMENT && field->components == 1,
                "TOPODENSITY field '", field_name, "' must be ELEMENT domain with 1 component");
            lc->density = field;
        });

        // Keyword-only command: no data segment needed
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
