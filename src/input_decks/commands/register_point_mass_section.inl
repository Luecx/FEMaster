// register_point_mass_section.inl — DSL registration for *POINTMASSSECTION

#include <array>
#include <string>

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_point_mass_section(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("POINTMASSSECTION", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign point-mass properties to an element set.");

        std::string elset;

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("ELSET")
                    .required()
                    .doc("Target element set")
        );

        command.on_enter([&](const fem::dsl::Keys& keys) {
            elset = keys.raw("ELSET");
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .allow_multiline()
                    .one<fem::Precision>().name("MASS").desc("Point mass")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 3>().name("INERTIA").desc("Rotary inertia Ix,Iy,Iz")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 3>().name("SPRING").desc("Translational spring constants")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .fixed<fem::Precision, 3>().name("ROTSPRING").desc("Rotational spring constants")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, &elset](fem::Precision mass,
                                       const std::array<fem::Precision, 3>& inertia_data,
                                       const std::array<fem::Precision, 3>& spring_data,
                                       const std::array<fem::Precision, 3>& rotary_data) {
                    fem::Vec3 inertia;
                    inertia << inertia_data[0], inertia_data[1], inertia_data[2];

                    fem::Vec3 spring;
                    spring << spring_data[0], spring_data[1], spring_data[2];

                    fem::Vec3 rotary;
                    rotary << rotary_data[0], rotary_data[1], rotary_data[2];

                    model.point_mass_section(elset, mass, inertia, spring, rotary);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
