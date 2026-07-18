// register_elastic.inl — DSL registration for *ELASTIC

#include <array>
#include <stdexcept>

#include "../../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../material/generalised_isotropic_elasticity.h"
#include "../../../material/isotropic_elasticity.h"
#include "../../../material/neo_hooke_elasticity.h"
#include "../../../material/orthotropic_elasticity.h"

namespace fem::io::reader::commands {

inline void register_elastic(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("ELASTIC", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc("Assign elastic properties to the active material.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE")
                    .optional("ISO")
                    .doc("Elasticity formulation (ISO/NEO_HOOKE/GENERALISED ISOTROPIC/ORTHO)")
                    .allowed({"ISO", "ISOTROPIC",
                              "NEOHOOKE", "NEO_HOOKE", "NEO-HOOKE", "NEOHOOKEAN", "NEO_HOOKEAN",
                              "GENERALISEDISOTROPIC", "GENERALISED_ISOTROPIC", "GENISO",
                              "ORTHO", "ORTHOTROPIC"})
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"ISO", "ISOTROPIC"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
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

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"NEOHOOKE", "NEO_HOOKE", "NEO-HOOKE", "NEOHOOKEAN", "NEO_HOOKEAN"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
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
                    material->set_elasticity<fem::material::NeoHookeElasticity>(E, nu);
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"GENERALISEDISOTROPIC", "GENERALISED_ISOTROPIC", "GENISO"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("E").desc("Young's modulus")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .one<fem::Precision>().name("NU").desc("Poisson ratio")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                    .one<fem::Precision>().name("G").desc("Independent shear modulus")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model](fem::Precision E, fem::Precision nu, fem::Precision G) {
                    auto material = model._data->materials.get();
                    if (!material) {
                        throw std::runtime_error("ELASTIC requires an active material context");
                    }
                    material->set_elasticity<fem::material::GeneralisedIsotropicElasticity>(E, nu, G);
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"ORTHO", "ORTHOTROPIC"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1).max(1))
                .pattern(fem::io::dsl::Pattern::make()
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
                    const fem::Precision nu31 = nu13 * E3 / E1;
                    material->set_elasticity<fem::material::OrthotropicElasticity>(
                        E1, E2, E3, G23, G13, G12, nu23, nu31, nu12);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
