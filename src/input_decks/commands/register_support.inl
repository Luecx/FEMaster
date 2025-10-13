// register_support.inl — registers *SUPPORT

#include <array>
#include <limits>
#include <stdexcept>
#include <string>

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_support(fem::dsl::Registry& registry, model::Model& model) {
    std::string orientation;

    registry.command("SUPPORT", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Define nodal supports via support collectors.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("SUPPORT_COLLECTOR").required().doc("Support collector name")
                .key("ORIENTATION").optional().doc("Optional orientation coordinate system")
        );

        command.on_enter([&](const fem::dsl::Keys& keys) {
            const std::string& collector = keys.raw("SUPPORT_COLLECTOR");
            orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            model._data->supp_cols.activate(collector);
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Node set or id")
                    .fixed<fem::Precision, 6>().name("DOF").desc("Support values for ux,uy,uz,rx,ry,rz")
                        .on_missing(std::numeric_limits<fem::Precision>::quiet_NaN())
                        .on_empty(std::numeric_limits<fem::Precision>::quiet_NaN())
                )
                .bind([&](const std::string& target,
                          const std::array<fem::Precision, 6>& values) {
                    fem::StaticVector<6> constraint;
                    for (int i = 0; i < 6; ++i) constraint(i) = values[i];

                    if (model._data->node_sets.has(target)) {
                        model.add_support(target, constraint, orientation);
                        return;
                    }

                    try {
                        const fem::ID id = static_cast<fem::ID>(std::stoi(target));
                        model.add_support(id, constraint, orientation);
                    } catch (const std::exception&) {
                        throw std::runtime_error("SUPPORT target '" + target + "' is neither a node set nor an id");
                    }
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
