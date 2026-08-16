// register_plastic.inl — DSL registration for *PLASTIC

#include <stdexcept>

#include "../../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../material/j2_plasticity.h"

namespace fem::io::reader::commands {

inline void register_plastic(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("PLASTIC", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("MATERIAL"));
        command.doc(
            "Assign small-strain J2 plasticity with tabulated isotropic hardening. "
            "Each data line contains yield stress and equivalent plastic strain; "
            "the first PEEQ must be zero."
        );

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE")
                    .optional("J2")
                    .doc("Plasticity formulation")
                    .allowed({"J2", "MISES", "VONMISES", "VON_MISES"})
        );

        command.on_enter([&model](const fem::io::dsl::Keys&) {
            auto material = model._data->materials.get();
            if (!material) {
                throw std::runtime_error("PLASTIC requires an active material context");
            }
            material->set_plasticity<fem::material::J2Plasticity>();
        });

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(1))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<fem::Precision>().name("YIELD_STRESS").desc("Current yield stress")
                    .one<fem::Precision>().name("PEEQ").desc("Equivalent plastic strain")
                )
                .bind([&model](fem::Precision yield_stress, fem::Precision peeq) {
                    auto material = model._data->materials.get();
                    if (!material || !material->has_plasticity()) {
                        throw std::runtime_error("PLASTIC requires an active plastic material context");
                    }
                    auto j2 = material->plasticity()->as<fem::material::J2Plasticity>();
                    if (!j2) {
                        throw std::runtime_error("PLASTIC data do not match active plasticity type");
                    }
                    j2->add_hardening_point(yield_stress, peeq);
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
