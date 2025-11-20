// register_beam_section.inl — DSL registration for *BEAMSECTION
#pragma once
/**
 * @file register_beam_section.inl
 * @brief Register *BEAMSECTION command.
 *
 * Assigns a beam section (material + profile) and a local section orientation vector.
 * The orientation is orthonormalized against the element axis internally; near-collinear
 * inputs trigger a robust fallback (implementation-defined).
 */

#include <array>
#include <memory>
#include <string>

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_beam_section(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("BEAMSECTION", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));

        // Keep this concise; no explicit “scope/keywords/data” lists here.
        command.doc(
            "Assign a beam section with a local orientation. "
            "The target element set must contain only beam elements. "
            "The provided direction n1 is normalized and orthogonalized internally; "
            "if n1 is (near-)collinear with the element axis, a stable fallback is applied."
        );

        // Persistent per-command state:
        auto material    = std::make_shared<std::string>();
        auto elset       = std::make_shared<std::string>();
        auto profile     = std::make_shared<std::string>();
        auto orientation = std::make_shared<fem::Vec3>(fem::Vec3::Zero());

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("MATERIAL").alternative("MAT").required()
                    .doc("Material name (must exist).")
                .key("ELSET").required()
                    .doc("Target element set; must contain only beam elements.")
                .key("PROFILE").required()
                    .doc("Section/profile identifier, e.g., IPE80, RECT_20x10, CHS60x3.")
        );

        command.on_enter([material, elset, profile, orientation](const fem::dsl::Keys& keys) {
            *material = keys.raw("MATERIAL");
            *elset    = keys.raw("ELSET");
            *profile  = keys.raw("PROFILE");
            orientation->setZero();
        });

        command.on_exit([&model, material, elset, profile, orientation](const fem::dsl::Keys&) {
            model.beam_section(*elset, *material, *profile, *orientation);
        });

        command.variant(
            fem::dsl::Variant::make()
                .doc("Optional data line with the section direction n1 = (N1_x, N1_y, N1_z).")
                .segment(
                    fem::dsl::Segment::make()
                        .range(fem::dsl::LineRange{}.min(0).max(1))
                        .pattern(
                            fem::dsl::Pattern::make()
                                .fixed<fem::Precision, 3>()
                                .name("N1")
                                .desc("Section orientation vector components: N1_x, N1_y, N1_z.")
                                .on_missing(fem::Precision{0})
                                .on_empty  (fem::Precision{0})
                        )
                        .bind([orientation](const std::array<fem::Precision, 3>& n1_data) {
                            fem::Vec3 n1; n1 << n1_data[0], n1_data[1], n1_data[2];
                            *orientation = n1;
                        })
                )
        );
    });
}

} // namespace fem::input_decks::commands
