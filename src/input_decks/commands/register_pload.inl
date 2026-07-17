// register_pload.inl — DSL registration for *PLOAD

#include <stdexcept>
#include <string>

#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_pload(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("PLOAD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Pressure loads on surfaces. AMPLITUDE may reference a time history that scales the pressure magnitude."
        );

        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR")
                    .doc("Target load collector that groups the loads")
                    .required()
                .key("AMPLITUDE").optional()
                    .doc("Optional time amplitude used to scale the pressure value")
        );

        command.on_enter([&model, amplitude](const fem::io::dsl::Keys& keys) {
            const std::string& collector = keys.raw("LOAD_COLLECTOR");
            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            model._data->load_cols.activate(collector);
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Surface set or surface id")
                    .one<fem::Precision>().name("P").desc("Pressure value")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, amplitude](const std::string& target, fem::Precision value) {
                    if (model._data->surface_sets.has(target)) {
                        model.add_pload(target, value, *amplitude);
                        return;
                    }

                    try {
                        const fem::ID id = static_cast<fem::ID>(std::stoi(target));
                        model.add_pload(id, value, *amplitude);
                    } catch (const std::exception&) {
                        throw std::runtime_error("PLOAD target '" + target + "' is neither a surface set nor an id");
                    }
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
