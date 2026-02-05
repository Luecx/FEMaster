// register_inertiaload.inl — DSL registration for *INERTIALOAD

#include <array>
#include <memory>
#include <stdexcept>
#include <string>

#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

inline void register_inertialload(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("INERTIALOAD", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Rigid-body inertial load on element sets. Provide CENTER, CENTER_ACC, OMEGA, ALPHA. "
            "TARGET must be an element set. Coriolis term is not included."
        );

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR")
                    .doc("Load collector that stores the inertial loads")
                    .required()
        );

        command.on_enter([&model](const fem::dsl::Keys& keys) {
            const std::string& collector = keys.raw("LOAD_COLLECTOR");
            model._data->load_cols.activate(collector);
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Element set name")
                    .fixed<fem::Precision, 3>().name("CENTER").desc("Center point x,y,z")
                    .fixed<fem::Precision, 3>().name("CENTER_ACC").desc("Center linear acceleration ax,ay,az")
                    .fixed<fem::Precision, 3>().name("OMEGA").desc("Angular velocity wx,wy,wz")
                    .fixed<fem::Precision, 3>().name("ALPHA").desc("Angular acceleration ax,ay,az")
                )
                .bind([&model](const std::string& target,
                               const std::array<fem::Precision, 3>& center,
                               const std::array<fem::Precision, 3>& center_acc,
                               const std::array<fem::Precision, 3>& omega,
                               const std::array<fem::Precision, 3>& alpha) {
                    fem::Vec3 c; c << center[0], center[1], center[2];
                    fem::Vec3 a0; a0 << center_acc[0], center_acc[1], center_acc[2];
                    fem::Vec3 w;  w  << omega[0], omega[1], omega[2];
                    fem::Vec3 al; al << alpha[0], alpha[1], alpha[2];

                    if (model._data->elem_sets.has(target)) {
                        model.add_inertialload(target, c, a0, w, al);
                        return;
                    }

                    throw std::runtime_error("INERTIALOAD target '" + target + "' must be an element set");
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands

