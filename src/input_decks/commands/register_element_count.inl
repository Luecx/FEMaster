// register_element_count.inl — first-stage ID counting for *ELEMENT

#include <array>
#include <functional>
#include <string>
#include <utility>

#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../dsl/registry.h"

namespace fem::input_decks::commands {

using ElementCountSink = std::function<void(fem::ID)>;

inline void register_element_count(fem::dsl::Registry& registry, ElementCountSink sink) {
    registry.command("ELEMENT", [sink = std::move(sink)](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Count element ids.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("ELSET")
                    .optional("EALL")
                    .doc("Element set to activate while reading")
                .key("TYPE")
                    .required()
                    .doc("Element topology/type")
                    .allowed({
                        "C3D4", "C3D5", "C3D6", "C3D8", "C3D10", "C3D15", "C3D20", "C3D20R",
                        "B33", "T3", "S3", "S4", "MITC4", "MITC4FRT", "S6", "S8", "QSPT"
                    })
        );

#define FEM_ADD_ELEMENT_COUNT_VARIANT(KEY, COUNT, DESC) \
        command.variant(fem::dsl::Variant::make() \
            .when(fem::dsl::Condition::key_equals("TYPE", {KEY})) \
            .segment(fem::dsl::Segment::make() \
                .range(fem::dsl::LineRange{}.min(1)) \
                .pattern(fem::dsl::Pattern::make() \
                    .allow_multiline() \
                    .one<fem::ID>().name("ID").desc("Element id") \
                    .fixed<fem::ID, COUNT>().name("N").desc(DESC) \
                ) \
                .bind([sink](fem::ID id, const std::array<fem::ID, COUNT>&) { \
                    sink(id); \
                }) \
            ) \
        );

        FEM_ADD_ELEMENT_COUNT_VARIANT("C3D4", 4, "C3D4 connectivity (4 nodes)");
        FEM_ADD_ELEMENT_COUNT_VARIANT("C3D6", 6, "C3D6 connectivity (6 nodes)");
        FEM_ADD_ELEMENT_COUNT_VARIANT("C3D8", 8, "C3D8 connectivity (8 nodes)");
        FEM_ADD_ELEMENT_COUNT_VARIANT("C3D10", 10, "C3D10 connectivity (10 nodes)");
        FEM_ADD_ELEMENT_COUNT_VARIANT("C3D15", 15, "C3D15 connectivity (15 nodes)");
        FEM_ADD_ELEMENT_COUNT_VARIANT("C3D20", 20, "C3D20 connectivity (20 nodes)");
        FEM_ADD_ELEMENT_COUNT_VARIANT("C3D20R", 20, "C3D20R connectivity (20 nodes)");
        FEM_ADD_ELEMENT_COUNT_VARIANT("T3", 2, "T3 connectivity (2 nodes)");
        FEM_ADD_ELEMENT_COUNT_VARIANT("S3", 3, "S3 connectivity (3 nodes)");
        FEM_ADD_ELEMENT_COUNT_VARIANT("S4", 4, "S4 connectivity (4 nodes)");
        FEM_ADD_ELEMENT_COUNT_VARIANT("QSPT", 4, "QSPT connectivity (4 nodes)");
        FEM_ADD_ELEMENT_COUNT_VARIANT("MITC4", 4, "MITC4 connectivity (4 nodes)");
        FEM_ADD_ELEMENT_COUNT_VARIANT("MITC4FRT", 4, "MITC4FRT connectivity (4 nodes)");
        FEM_ADD_ELEMENT_COUNT_VARIANT("S6", 6, "S6 connectivity (6 nodes)");
        FEM_ADD_ELEMENT_COUNT_VARIANT("S8", 8, "S8 connectivity (8 nodes)");

#undef FEM_ADD_ELEMENT_COUNT_VARIANT

        command.variant(fem::dsl::Variant::make()
            .when(fem::dsl::Condition::key_equals("TYPE", {"B33"}))
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<fem::ID>().name("ID").desc("Element id")
                    .fixed<fem::ID, 3>()
                        .name("N")
                        .desc("B33 connectivity (2 nodes plus optional orientation node)")
                        .on_missing(fem::ID{-1})
                )
                .bind([sink](fem::ID id, const std::array<fem::ID, 3>&) {
                    sink(id);
                })
            )
        );

        command.variant(fem::dsl::Variant::make()
            .when(fem::dsl::Condition::key_equals("TYPE", {"C3D5"}))
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .allow_multiline()
                    .one<fem::ID>().name("ID").desc("Element id")
                    .fixed<fem::ID, 5>().name("N").desc("C3D5 connectivity (5 nodes)")
                )
                .bind([sink](fem::ID id, const std::array<fem::ID, 5>&) {
                    sink(id);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
