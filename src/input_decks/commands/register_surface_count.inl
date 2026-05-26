// register_surface_count.inl — first-stage ID counting for *SURFACE

#include <functional>
#include <string>
#include <utility>

#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"

namespace fem::input_decks::commands {

using SurfaceCountSink = std::function<void(fem::ID)>;

inline void register_surface_count(fem::dsl::Registry& registry, SurfaceCountSink sink) {
    registry.command("SURFACE", [sink = std::move(sink)](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Count explicit surface ids.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("SFSET").alternative("NAME").optional("SFALL")
                    .doc("Set name to activate/create (default: SFALL).")
        );

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<fem::ID>().name("ID").desc("Surface id")
                    .one<fem::ID>().name("ELEM_ID").desc("Element id")
                    .one<std::string>().name("SIDE").desc("Face id, e.g. S1 or 1")
                )
                .bind([sink](fem::ID id, fem::ID, const std::string&) {
                    sink(id);
                })
            )
        );

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("ELSET name or element id")
                    .one<std::string>().name("SIDE").desc("Face id, e.g. S1 or 1")
                )
            )
        );
    });
}

} // namespace fem::input_decks::commands
