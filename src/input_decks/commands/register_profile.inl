// ===== FILE: ./include/fem/input_decks/commands/register_profile.inl =====
#include <memory>
#include <string>

#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_profile(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("PROFILE", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Define beam profile properties in this order: A, Iy, Iz, Jt, Iyz, ey, ez, refy, refz. "
                    "Only the first 4 are required. "
                    "Convention: Iyz = integral_A(y*z*dA), i.e. without a leading minus sign.");

        // Persistent state for this command's handlers:
        auto profile_name = std::make_shared<std::string>();

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("PROFILE")
                    .alternative("NAME")
                    .required()
                    .doc("Identifier of the profile")
        );

        // Capture by value (shared_ptr) so it stays valid
        command.on_enter([profile_name](const fem::dsl::Keys& keys) {
            *profile_name = keys.raw("PROFILE");  // spec canonicalizes NAME→PROFILE
        });

        // Single variant with optional trailing values.
        // Order: A, Iy, Iz, Jt, Iyz, ey, ez, refy, refz
        // Missing/empty trailing values default to 0.
        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .fixed<fem::Precision, 1>().name("A").desc("Cross-section area A")
                    .fixed<fem::Precision, 1>().name("IY").desc("Second moment of area about local y-axis (Iy)")
                    .fixed<fem::Precision, 1>().name("IZ").desc("Second moment of area about local z-axis (Iz)")
                    .fixed<fem::Precision, 1>().name("JT").desc("Torsional constant (Jt)")
                    .fixed<fem::Precision, 1>().name("IYZ").desc("Product of inertia: Iyz = integral_A(y*z*dA), no minus sign")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 1>().name("EY").desc("Offset in local y: ey = y(SP) - y(SMP)")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 1>().name("EZ").desc("Offset in local z: ez = z(SP) - z(SMP)")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 1>().name("REFY").desc("Reference-line offset in local y: refy = y(REF) - y(SMP)")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 1>().name("REFZ").desc("Reference-line offset in local z: refz = z(REF) - z(SMP)")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, profile_name](fem::Precision A,
                                             fem::Precision IY,
                                             fem::Precision IZ,
                                             fem::Precision JT,
                                             fem::Precision IYZ,
                                             fem::Precision EY,
                                             fem::Precision EZ,
                                             fem::Precision REFY,
                                             fem::Precision REFZ) {
                    model._data->profiles.activate(*profile_name,
                                                   A, IY, IZ, JT, IYZ,
                                                   EY, EZ, REFY, REFZ);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
