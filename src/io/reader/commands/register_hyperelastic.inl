// register_hyperelastic.inl — registers *HYPERELASTIC

#include <stdexcept>

#include "../../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../material/neo_hooke_elasticity.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

inline void register_hyperelastic(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("HYPERELASTIC", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign an Abaqus-style hyperelastic model to the active material.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .flag("NEOHOOKE")
                    .alternative("NEO_HOOKE")
                    .alternative("NEO-HOOKE")
                    .required()
                    .doc("Use the Neo-Hooke potential with C10 and D1 parameters.")
        );

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("C10").desc("Neo-Hooke deviatoric coefficient")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .one<fem::Precision>().name("D1").desc("Neo-Hooke compressibility coefficient")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model](fem::Precision c10, fem::Precision d1) {
                    auto material = model._data->materials.get();
                    if (!material) {
                        throw std::runtime_error("HYPERELASTIC requires an active material context");
                    }
                    material->set_elasticity<fem::material::NeoHookeElasticity>(c10, d1);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
