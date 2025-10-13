// register_elastic.inl — DSL registration for *ELASTIC

#include <array>
#include <stdexcept>

#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../material/isotropic_elasticity.h"
#include "../../material/orthotropic_elasticity.h"

namespace fem::input_decks::commands {

inline void register_elastic(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("ELASTIC", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign elastic properties to the active material.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("TYPE")
                    .required()
                    .doc("Elasticity formulation (ISO/ORTHO)")
                    .allowed({"ISO", "ISOTROPIC", "ORTHO", "ORTHOTROPIC"})
        );

        command.variant(fem::dsl::Variant::make()
            .when(fem::dsl::Condition::key_equals("TYPE", {"ISO", "ISOTROPIC"}))
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<fem::Precision>().name("E").desc("Young's modulus")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .one<fem::Precision>().name("NU").desc("Poisson ratio")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model](fem::Precision E, fem::Precision nu) {
                    auto material = model._data->materials.get();
                    if (!material) {
                        throw std::runtime_error("ELASTIC requires an active material context");
                    }
                    material->set_elasticity<fem::material::IsotropicElasticity>(E, nu);
                })
            )
        );

        command.variant(fem::dsl::Variant::make()
            .when(fem::dsl::Condition::key_equals("TYPE", {"ORTHO", "ORTHOTROPIC"}))
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<fem::Precision>().name("E1").desc("Young's modulus in 1-direction")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .one<fem::Precision>().name("E2").desc("Young's modulus in 2-direction")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .one<fem::Precision>().name("E3").desc("Young's modulus in 3-direction")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .one<fem::Precision>().name("G23").desc("Shear modulus G23")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .one<fem::Precision>().name("G13").desc("Shear modulus G13")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .one<fem::Precision>().name("G12").desc("Shear modulus G12")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .one<fem::Precision>().name("NU23").desc("Poisson ratio nu23")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .one<fem::Precision>().name("NU13").desc("Poisson ratio nu13")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .one<fem::Precision>().name("NU12").desc("Poisson ratio nu12")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model](fem::Precision E1,
                               fem::Precision E2,
                               fem::Precision E3,
                               fem::Precision G23,
                               fem::Precision G13,
                               fem::Precision G12,
                               fem::Precision nu23,
                               fem::Precision nu13,
                               fem::Precision nu12) {
                    auto material = model._data->materials.get();
                    if (!material) {
                        throw std::runtime_error("ELASTIC requires an active material context");
                    }
                    material->set_elasticity<fem::material::OrthotropicElasticity>(
                        E1, E2, E3, G23, G13, G12, nu23, nu13, nu12);
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
