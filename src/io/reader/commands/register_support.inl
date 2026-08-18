// register_support.inl — registers *SUPPORT

#include <array>
#include <charconv>
#include <limits>
#include <memory>
#include <string>
#include <system_error>

#include "../../../core/logging.h"
#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

inline void register_support(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("SUPPORT", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Define nodal supports via support collectors.");

        // Persistent per-command state (captured by value in lambdas)
        auto orientation = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("SUPPORT_COLLECTOR").required().doc("Support collector name")
                .key("ORIENTATION").optional().doc("Optional orientation coordinate system")
        );

        command.on_enter([orientation, &model](const fem::io::dsl::Keys& keys) {
            // NOTE: raw() presumably returns by value; if not, copy explicitly.
            const std::string collector = keys.raw("SUPPORT_COLLECTOR");
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            model._data->supp_cols.activate(collector);
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Node set or id")
                    .fixed<fem::Precision, 6>().name("DOF").desc("Support values for ux,uy,uz,rx,ry,rz")
                        .on_missing(std::numeric_limits<fem::Precision>::quiet_NaN())
                        .on_empty  (std::numeric_limits<fem::Precision>::quiet_NaN())
                )
                .bind([&model, orientation](const std::string& target,
                                            const std::array<fem::Precision, 6>& values) {
                    fem::StaticVector<6> constraint;
                    for (int i = 0; i < 6; ++i) constraint(i) = values[i];

                    if (model._data->node_sets.has(target)) {
                        model.add_support(target, constraint, *orientation);
                        return;
                    }

                    fem::ID id{};
                    const char* begin = target.data();
                    const char* end   = begin + target.size();
                    const auto [ptr, ec] = std::from_chars(begin, end, id);
                    logging::error(ec == std::errc{} && ptr == end,
                        "SUPPORT target '", target, "' is neither a node set nor an id");
                    model.add_support(id, constraint, *orientation);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
