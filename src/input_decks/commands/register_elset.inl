// register_elset.inl — DSL registration for *ELSET

#include <array>
#include <limits>
#include <string>

#include "../../core/types_num.h"
#include "../../dsl/keyword.h"

namespace fem::input_decks::commands {

inline void register_elset(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("ELSET", [&](fem::dsl::Command& command) {
        command.doc("Define named element sets.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("ELSET")
                    .alternative("NAME")
                    .required()
                    .doc("Name of the element set to activate or create")
                .flag("GENERATE")
                    .doc("Interpret rows as start,end[,increment] ranges")
        );

        command.on_enter([&model](const fem::dsl::Keys& keys) {
            const std::string& name = keys.raw("ELSET");
            model._data->elem_sets.activate(name);
        });

        const fem::ID missing_id = std::numeric_limits<fem::ID>::min();

        command.variant(fem::dsl::Variant::make()
            .when(fem::dsl::Condition::key_present("GENERATE"))
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<fem::ID>().name("START").desc("First element id")
                    .one<fem::ID>().name("END").desc("Last element id")
                    .one<fem::ID>().name("INC").desc("Increment").on_missing(fem::ID{1}).on_empty(fem::ID{1})
                )
                .bind([&model](fem::ID first, fem::ID last, fem::ID inc) {
                    if (inc == 0) {
                        inc = 1;
                    }

                    auto set = model._data->elem_sets.get();
                    if (!set) {
                        return;
                    }

                    if ((inc > 0 && first > last) || (inc < 0 && first < last)) {
                        return;
                    }

                    for (fem::ID id = first;; id += inc) {
                        set->add(id);

                        const fem::ID next = static_cast<fem::ID>(id + inc);
                        const bool stop = (inc > 0) ? next > last : next < last;
                        if (stop) {
                            break;
                        }
                    }
                })
            )
        );

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(0))
                .pattern(fem::dsl::Pattern::make()
                    .fixed<fem::ID, 32>().name("ID").desc("Element ids (up to 32 per line)")
                        .on_missing(missing_id)
                        .on_empty(missing_id)
                )
                .bind([&model, missing_id](const std::array<fem::ID, 32>& ids) {
                    auto set = model._data->elem_sets.get();
                    if (!set) {
                        return;
                    }
                    for (fem::ID id : ids) {
                        if (id == missing_id) {
                            continue;
                        }
                        set->add(id);
                    }
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
