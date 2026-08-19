/**
 * @file register_profile.inl
 * @brief Registers object-based beam profile definitions for FEMaster decks.
 *
 * Beam cross-section constants are parsed into a complete `Profile` object whose
 * name is part of the object itself. The parser registers that object through
 * `Model::add_profile()` instead of asking the model facade to construct a
 * profile from duplicated name and scalar arguments.
 *
 * The product-of-inertia convention remains
 * `Iyz = integral_A(y*z*dA)` without a leading minus sign.
 *
 * @see Profile
 * @see model::Model::add_profile
 *
 * @author Finn Eggers
 * @date 18.08.2026
 */

#include <memory>
#include <string>

#include "../../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"
#include "../../../section/profile.h"

namespace fem::io::reader::commands {

inline void register_profile(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("PROFILE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Define beam profile properties in this order: A, Iy, Iz, Jt, Iyz, ey, ez, refy, refz. "
                    "Only the first 4 are required. "
                    "Convention: Iyz = integral_A(y*z*dA), i.e. without a leading minus sign.");

        auto profile_name = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("PROFILE")
                    .alternative("NAME")
                    .required()
                    .doc("Identifier of the profile")
        );

        command.on_enter([profile_name](const fem::io::dsl::Keys& keys) {
            *profile_name = keys.raw("PROFILE");
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
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
                .bind([&model, profile_name](fem::Precision area,
                                             fem::Precision inertia_y,
                                             fem::Precision inertia_z,
                                             fem::Precision torsion,
                                             fem::Precision product_yz,
                                             fem::Precision offset_y,
                                             fem::Precision offset_z,
                                             fem::Precision reference_y,
                                             fem::Precision reference_z) {
                    model.add_profile(std::make_shared<Profile>(
                        *profile_name,
                        area,
                        inertia_y,
                        inertia_z,
                        torsion,
                        product_yz,
                        offset_y,
                        offset_z,
                        reference_y,
                        reference_z
                    ));
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands