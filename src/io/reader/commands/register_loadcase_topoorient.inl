/**
 * @file register_loadcase_topoorient.inl
 * @brief Registers element orientations for topology-weighted statics.
 *
 * `TOPOORIENT` resolves a named three-component element field and assigns it to
 * the active `LinearStaticTopo` analysis. The field supplies the element-wise
 * orientation data required by orientation-dependent topology calculations.
 *
 * Domain and component-count checks are performed while parsing. Interpretation
 * of the stored orientation vectors remains within the load-case formulation.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "../parser.h"

#include "../../../core/logging.h"
#include "../../../loadcase/linear_static_topo.h"

namespace fem::io::reader::commands {

inline void register_loadcase_topoorient(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("TOPOORIENT", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Select the orientation field for LINEARSTATICTOPO loadcases.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("FIELD").required().doc("Orientation field name (ELEMENT, 3 components)")
                .alternative("NAME")
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto* lc = dynamic_cast<loadcase::LinearStaticTopo*>(parser.active_loadcase());
            logging::error(lc != nullptr,
                "TOPOORIENT only valid for LINEARSTATICTOPO loadcases");

            const std::string field_name = keys.raw("FIELD");
            auto field = parser.model()._data->get_field(field_name);
            logging::error(field != nullptr,
                "TOPOORIENT field '", field_name, "' does not exist");
            logging::error(field->domain == model::FieldDomain::ELEMENT && field->components == 3,
                "TOPOORIENT field '", field_name, "' must be ELEMENT domain with 3 components");
            lc->orientation = field;
        });

        // Keyword-only command: no data segment needed
        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
