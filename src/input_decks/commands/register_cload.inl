// register_cload.inl — DSL registration for *CLOAD

#include <array>
#include <stdexcept>
#include <string>

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_cload(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("CLOAD", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Concentrated nodal loads (forces and moments).");

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
                    .one<std::string>().name("TARGET").desc("Node set or node id")
                    .fixed<fem::Precision, 6>().name("LOAD").desc("Fx,Fy,Fz,Mx,My,Mz")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model](const std::string& target, const std::array<fem::Precision, 6>& values) {
                    fem::Vec6 load;
                    load << values[0], values[1], values[2], values[3], values[4], values[5];

                    if (model._data->node_sets.has(target)) {
                        model.add_cload(target, load);
                        return;
                    }

                    try {
                        const fem::ID id = static_cast<fem::ID>(std::stoi(target));
                        model.add_cload(id, load);
                    } catch (const std::exception&) {
                        throw std::runtime_error("CLOAD target '" + target + "' is neither a node set nor an id");
                    }
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
