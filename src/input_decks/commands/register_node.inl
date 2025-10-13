// register_node.inl — DSL registration for *NODE

#include <string>

#include "../../core/types_num.h"

#include "../../dsl/keyword.h"

namespace fem::input_decks::commands {

inline void register_node(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("NODE", [&](fem::dsl::Command& command) {
        command.doc("Define nodes with optional coordinates and assign them to a node set.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("NSET")
                    .optional("NALL")
                    .doc("Target node set to populate (defaults to NALL)")
        );

        command.on_enter([&model](const fem::dsl::Keys& keys) {
            const std::string& target = keys.raw("NSET");
            model._data->node_sets.activate(target);
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(0))
                .pattern(fem::dsl::Pattern::make()
                    .one<int>().name("ID").desc("Node identifier")
                    .one<fem::Precision>().name("X").desc("X-coordinate").on_empty(fem::Precision{0}).on_missing(fem::Precision{0})
                    .one<fem::Precision>().name("Y").desc("Y-coordinate").on_empty(fem::Precision{0}).on_missing(fem::Precision{0})
                    .one<fem::Precision>().name("Z").desc("Z-coordinate").on_empty(fem::Precision{0}).on_missing(fem::Precision{0})
                )
                .bind([&model](int id, fem::Precision x, fem::Precision y, fem::Precision z) {
                    model.set_node(id, x, y, z);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
