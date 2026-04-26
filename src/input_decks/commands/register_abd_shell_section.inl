// register_abd_shell_section.inl - DSL registration for *ABDSHELLSECTION

#include <array>
#include <memory>
#include <string>

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_abd_shell_section(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("ABDSHELLSECTION", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign a shell section from explicit ABD (6x6) and transverse shear (2x2) matrices.");

        auto elset = std::make_shared<std::string>();
        auto orientation = std::make_shared<std::string>();
        auto thickness = std::make_shared<fem::Precision>(fem::Precision{1});

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("ELSET").required().doc("Target element set")
                .key("ORIENTATION").optional().doc("Optional coordinate system for shell n1/n2 material/resultant axes")
                .key("THICKNESS").optional().doc("Shell thickness used for geometric volume/load integration")
        );

        command.on_enter([elset, orientation, thickness](const fem::dsl::Keys& keys) {
            *elset = keys.raw("ELSET");
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            *thickness = keys.has("THICKNESS") ? keys.get<fem::Precision>("THICKNESS") : fem::Precision{1};
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(8))
                .pattern(fem::dsl::Pattern::make()
                    .allow_multiline()
                    .fixed<fem::Precision, 40>().name("DATA").desc("ABD row-major (36 values), then shear row-major (4 values)")
                )
                .bind([&model, elset, orientation, thickness](const std::array<fem::Precision, 40>& vals) {
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

                    model.abd_shell_section(*elset, *thickness, abd, shear, *orientation);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
