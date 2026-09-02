/**
 * @file register_plastic.cpp
 * @brief Registers tabulated isotropic J2 plasticity for material definitions.
 *
 * `PLASTIC` augments an already defined isotropic `ELASTIC` law. During semantic
 * execution the existing `IsotropicElasticity` is replaced by
 * `IsotropicJ2Elasticity` while preserving Young's modulus and Poisson's ratio.
 * Every data line then appends one `(yield stress, equivalent plastic strain)`
 * hardening point to the J2 model.
 *
 * Because material commands are executed in explicit semantic order, `ELASTIC`
 * is processed before `PLASTIC` even when the keywords appear in the opposite
 * order in the source deck.
 *
 * @see material::IsotropicElasticity
 * @see material::IsotropicJ2Elasticity
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../../../core/logging.h"
#include "../../../material/isotropic_elasticity.h"
#include "../../../material/isotropic_j2_elasticity.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

/**
 * Registers isotropic J2 plasticity for the active material.
 *
 * The command intentionally accepts only a preceding `IsotropicElasticity`.
 * This prevents arbitrary combinations of unrelated elastic and plastic laws.
 * Data follow the Abaqus-style order `yield stress, equivalent plastic strain`.
 */
void register_plastic(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("PLASTIC", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc(
            "Replace isotropic ELASTIC by isotropic J2 elastoplasticity with "
            "tabulated isotropic hardening."
        );

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE")
                    .optional("J2")
                    .allowed({"J2", "MISES", "VONMISES", "VON_MISES"})
                    .doc("Plasticity formulation")
        );

        command.on_enter([&model](const fem::io::dsl::Keys&) {
            auto material = model._data->materials.get();
            logging::error(material != nullptr,
                           "PLASTIC requires an active material context");
            logging::error(material->has_elasticity(),
                           "PLASTIC requires a preceding isotropic ELASTIC definition");

            const auto elasticity = material->elasticity();
            const auto* isotropic = elasticity->as<fem::material::IsotropicElasticity>();
            logging::error(isotropic != nullptr,
                           "PLASTIC currently supports only isotropic ELASTIC material definitions");

            const fem::Precision youngs  = isotropic->youngs;
            const fem::Precision poisson = isotropic->poisson;
            material->set_elasticity<fem::material::IsotropicJ2Elasticity>(youngs, poisson);
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("YIELD_STRESS").desc("Uniaxial yield stress")
                    .one<fem::Precision>().name("PEEQ").desc("Equivalent plastic strain")
                )
                .bind([&model](fem::Precision yield_stress, fem::Precision peeq) {
                    auto material = model._data->materials.get();
                    logging::error(material != nullptr && material->has_elasticity(),
                                   "PLASTIC data require an active material context");

                    auto* j2 = material->elasticity()->as<fem::material::IsotropicJ2Elasticity>();
                    logging::error(j2 != nullptr,
                                   "PLASTIC data do not match the active material law");
                    j2->add_yield_point(yield_stress, peeq);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
