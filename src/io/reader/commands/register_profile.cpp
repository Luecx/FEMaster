/**
 * @file register_profile.cpp
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

#include "register_functions.h"
#include "../../dsl/registry.h"

#include <memory>
#include <string>

#include "../../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"
#include "../../../section/profile.h"

namespace fem::io::reader::commands {

void register_profile(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("PROFILE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Define beam profile properties in this order: A, Iy, Iz, Jt, Iyz, ey, ez, refy, refz, Asy, Asz. "
                    "Only the first 4 are required; omitted shear areas default to 5A/6. "
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

        command.data(
            fem::io::dsl::Pattern::make()
                .fixed<fem::Precision, 1>("A", "Cross-section area A")
                .fixed<fem::Precision, 1>("IY", "Second moment of area about local y-axis (Iy)")
                .fixed<fem::Precision, 1>("IZ", "Second moment of area about local z-axis (Iz)")
                .fixed<fem::Precision, 1>("JT", "Torsional constant (Jt)")
                .fixed<fem::Precision, 1>("IYZ", "Product of inertia: Iyz = integral_A(y*z*dA), no minus sign")
                    .defaults(fem::Precision{0})
                .fixed<fem::Precision, 1>("EY", "Offset in local y: ey = y(SP) - y(SMP)")
                    .defaults(fem::Precision{0})
                .fixed<fem::Precision, 1>("EZ", "Offset in local z: ez = z(SP) - z(SMP)")
                    .defaults(fem::Precision{0})
                .fixed<fem::Precision, 1>("REFY", "Reference-line offset in local y: refy = y(REF) - y(SMP)")
                    .defaults(fem::Precision{0})
                .fixed<fem::Precision, 1>("REFZ", "Reference-line offset in local z: refz = z(REF) - z(SMP)")
                    .defaults(fem::Precision{0})
                .fixed<fem::Precision, 1>("ASY", "Effective transverse shear area in local y (default 5A/6)")
                    .defaults(fem::Precision{0})
                .fixed<fem::Precision, 1>("ASZ", "Effective transverse shear area in local z (default 5A/6)")
                    .defaults(fem::Precision{0}),
            [&model, profile_name](fem::Precision area,
                                         fem::Precision inertia_y,
                                         fem::Precision inertia_z,
                                         fem::Precision torsion,
                                         fem::Precision product_yz,
                                         fem::Precision offset_y,
                                         fem::Precision offset_z,
                                         fem::Precision reference_y,
                                         fem::Precision reference_z,
                                         fem::Precision shear_area_y,
                                         fem::Precision shear_area_z) {
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
                    reference_z,
                    shear_area_y,
                    shear_area_z
                ));
            },
            fem::io::dsl::LineRange{}.min(1).max(1)
        );
    });
}

} // namespace fem::io::reader::commands
