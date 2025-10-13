// register_vload.inl — DSL registration for *VLOAD

#include <array>
#include <stdexcept>
#include <string>

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_vload(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("VLOAD", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Volume loads applied to element sets or individual elements.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR")
                    .doc("Target load collector that groups the loads")
                    .required()
        );

        command.on_enter([&model](const fem::dsl::Keys& keys) {
            const std::string& collector = keys.raw("LOAD_COLLECTOR");
            model._data->load_cols.activate(collector);
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Element set or element id")
                    .fixed<fem::Precision, 3>().name("LOAD").desc("Fx,Fy,Fz")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model](const std::string& target, const std::array<fem::Precision, 3>& values) {
                    fem::Vec3 load;
                    load << values[0], values[1], values[2];

                    if (model._data->elem_sets.has(target)) {
                        model.add_vload(target, load);
                        return;
                    }

                    try {
                        const fem::ID id = static_cast<fem::ID>(std::stoi(target));
                        model.add_vload(id, load);
                    } catch (const std::exception&) {
                        throw std::runtime_error("VLOAD target '" + target + "' is neither an element set nor an id");
                    }
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
