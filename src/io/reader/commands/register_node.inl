// register_node.inl — DSL registration for *NODE

#include <string>

#include "../../../core/types_num.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

inline void register_node(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("NODE", [&](fem::io::dsl::Command& command) {
        command.doc("Define nodes with optional coordinates and assign them to a node set.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NSET")
                    .optional("NALL")
                    .doc("Target node set to populate (defaults to NALL)")
        );

        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            model._data->node_sets.activate(keys.raw("NSET"));
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::ID>().name("ID").desc("Node identifier")
                    .one<fem::Precision>().name("X").desc("X-coordinate").on_empty(fem::Precision{0}).on_missing(fem::Precision{0})
                    .one<fem::Precision>().name("Y").desc("Y-coordinate").on_empty(fem::Precision{0}).on_missing(fem::Precision{0})
                    .one<fem::Precision>().name("Z").desc("Z-coordinate").on_empty(fem::Precision{0}).on_missing(fem::Precision{0})
                )
                .bind([&model](fem::ID id, fem::Precision x, fem::Precision y, fem::Precision z) {
                    model.set_node(id, x, y, z);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
