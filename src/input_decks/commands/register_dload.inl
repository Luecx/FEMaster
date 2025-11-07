// register_dload.inl — DSL registration for *DLOAD

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

inline void register_dload(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("DLOAD", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc(
            "Apply distributed surface tractions to surfaces or surface sets. Each line provides a target followed by "
            "three components (Fx,Fy,Fz). When ORIENTATION is specified the components are measured in that local frame; "
            "otherwise they are interpreted in the global Cartesian basis. AMPLITUDE references a time history that scales "
            "the traction. Contributions from multiple records add linearly in the active load collector."
        );

        auto orientation = std::make_shared<std::string>();
        auto amplitude = std::make_shared<std::string>();

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("LOAD_COLLECTOR")
                    .doc("Name of the load collector that aggregates the distributed loads")
                    .required()
                .key("ORIENTATION").optional()
                    .doc("Optional coordinate system describing the local traction directions")
                .key("AMPLITUDE").optional()
                    .doc("Optional time amplitude that scales the traction vector")
        );

        command.on_enter([orientation, amplitude, &model](const fem::dsl::Keys& keys) {
            const std::string& collector = keys.raw("LOAD_COLLECTOR");
            *orientation = keys.has("ORIENTATION") ? keys.raw("ORIENTATION") : std::string{};
            *amplitude = keys.has("AMPLITUDE") ? keys.raw("AMPLITUDE") : std::string{};
            model._data->load_cols.activate(collector);
        });

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(1))
                .pattern(fem::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("Surface set or surface id")
                    .fixed<fem::Precision, 3>().name("LOAD").desc("Fx,Fy,Fz components, local if ORIENTATION is set")
                        .on_missing(fem::Precision{0}).on_empty(fem::Precision{0})
                )
                .bind([&model, orientation, amplitude](const std::string& target, const std::array<fem::Precision, 3>& values) {
                    fem::Vec3 load;
                    load << values[0], values[1], values[2];

                    if (model._data->surface_sets.has(target)) {
                        model.add_dload(target, load, *orientation, *amplitude);
                        return;
                    }

                    try {
                        const fem::ID id = static_cast<fem::ID>(std::stoi(target));
                        model.add_dload(id, load, *orientation, *amplitude);
                    } catch (const std::exception&) {
                        throw std::runtime_error("DLOAD target '" + target + "' is neither a surface set nor an id");
                    }
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
