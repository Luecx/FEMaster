// register_node_count.inl — first-stage ID counting for *NODE

#include <algorithm>
#include <functional>
#include <utility>

#include "../../core/types_num.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"

namespace fem::input_decks::commands {

using NodeCountSink = std::function<void(fem::ID)>;

inline void register_node_count(fem::dsl::Registry& registry, NodeCountSink sink) {
    registry.command("NODE", [sink = std::move(sink)](fem::dsl::Command& command) {
        command.doc("Count node ids.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("NSET")
                    .optional("NALL")
                    .doc("Target node set to populate (defaults to NALL)")
        );

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(0))
                .pattern(fem::dsl::Pattern::make()
                    .one<fem::ID>().name("ID").desc("Node identifier")
                    .one<fem::Precision>().name("X").desc("X-coordinate").on_empty(fem::Precision{0}).on_missing(fem::Precision{0})
                    .one<fem::Precision>().name("Y").desc("Y-coordinate").on_empty(fem::Precision{0}).on_missing(fem::Precision{0})
                    .one<fem::Precision>().name("Z").desc("Z-coordinate").on_empty(fem::Precision{0}).on_missing(fem::Precision{0})
                )
                .bind([sink](fem::ID id, fem::Precision, fem::Precision, fem::Precision) {
                    sink(id);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
