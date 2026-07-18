// register_point_mass.inl — DSL registration for *POINTMASS (feature on NSET)

#include <array>
#include <memory>
#include <string>

#include "../../../core/types_eig.h"
#include "../../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

inline void register_point_mass(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("POINTMASS", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Assign point-mass feature to a node set.");

        auto nset = std::make_shared<std::string>();

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NSET").required().doc("Target node set")
        );

        command.on_enter([nset](const fem::io::dsl::Keys& keys) {
            *nset = keys.raw("NSET");
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
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
                .bind([&model, nset](
                          fem::Precision mass,
                          const std::array<fem::Precision, 3>& inertia_data,
                          const std::array<fem::Precision, 3>& spring_data,
                          const std::array<fem::Precision, 3>& rotary_data) {
                    fem::Vec3 inertia{inertia_data[0], inertia_data[1], inertia_data[2]};
                    fem::Vec3 spring {spring_data[0],  spring_data[1],  spring_data[2]};
                    fem::Vec3 rotary {rotary_data[0],  rotary_data[1],  rotary_data[2]};
                    model.add_point_mass_feature(*nset, mass, inertia, spring, rotary);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands

