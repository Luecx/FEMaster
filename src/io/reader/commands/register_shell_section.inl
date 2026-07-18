// register_shell_section.inl — DSL registration for *SHELLSECTION

#include <memory>
#include <array>
#include <string>

#include "../../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

inline void register_shell_section(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("SHELLSECTION", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign a shell section to an element set.");

        // Persistent per-command state
        auto material    = std::make_shared<std::string>();
        auto elset       = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();
        auto type        = std::make_shared<std::string>();
        auto thickness   = std::make_shared<fem::Precision>(fem::Precision(1));

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE").optional("INTEGRATED").allowed({"INTEGRATED", "ABD"}).doc("Section formulation")
                .key("MATERIAL").alternative("MAT").optional().doc("Material name")
                .key("ELSET").required().doc("Target element set")
                .key("THICKNESS").optional("1.0").doc("Shell thickness used for mass/geometric thickness")
                .key("ORIENTATION").optional().doc("Optional coordinate system for shell n1/n2 material/resultant axes")
        );

        // Capture shared_ptrs BY VALUE so they're valid later
        command.on_enter([material, elset, orientation, type, thickness](const fem::io::dsl::Keys& keys) {
            *type        = keys.raw("TYPE");
            *material    = keys.has("MATERIAL") ? keys.raw("MATERIAL") : std::string{};
            *elset       = keys.raw("ELSET");
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            *thickness   = keys.get<fem::Precision>("THICKNESS");
        });

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"INTEGRATED"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("THICKNESS").desc("Shell thickness")
                        .on_missing(fem::Precision{1}).on_empty(fem::Precision{1})
                )
                .bind([&model, material, elset, orientation](fem::Precision thickness) {
                    if (material->empty()) {
                        throw std::runtime_error("SHELLSECTION TYPE=INTEGRATED requires MATERIAL");
                    }
                    model.shell_section(*elset, *material, thickness, *orientation);
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"ABD"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(8))
                .pattern(fem::io::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 40>().name("DATA").desc("ABD row-major (36 values), then shear row-major (4 values)")
                )
                .bind([&model, material, elset, orientation, thickness](const std::array<fem::Precision, 40>& vals) {
                    StaticMatrix<6, 6> abd;
                    StaticMatrix<2, 2> shear;
                    for (Index i = 0; i < 6; ++i) {
                        for (Index j = 0; j < 6; ++j) {
                            abd(i, j) = vals[6 * i + j];
                        }
                    }
                    for (Index i = 0; i < 2; ++i) {
                        for (Index j = 0; j < 2; ++j) {
                            shear(i, j) = vals[36 + 2 * i + j];
                        }
                    }

                    model.shell_section_abd(*elset, *material, *thickness, abd, shear, *orientation);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
