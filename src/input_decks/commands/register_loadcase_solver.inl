// register_loadcase_solver.inl — registers SOLVER within *LOADCASE

#include <stdexcept>
#include <string>

#include "../parser.h"

#include "../../loadcase/linear_buckling.h"
#include "../../loadcase/linear_static.h"
#include "../../loadcase/linear_static_topo.h"

namespace fem::input_decks::commands {

inline void register_loadcase_solver(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("SOLVER", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("LOADCASE"));
        command.doc("Configure solver options for the active loadcase.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("DEVICE").optional("CPU").allowed({"CPU", "GPU"})
                .key("METHOD").optional("DIRECT").allowed({"DIRECT", "INDIRECT"})
        );

        command.on_enter([&parser](const fem::dsl::Keys& keys) {
            auto* base = parser.active_loadcase();
            if (!base) {
                throw std::runtime_error("SOLVER must appear inside *LOADCASE");
            }

            const auto device = keys.raw("DEVICE");
            const auto method = keys.raw("METHOD");

            auto configure = [&](auto* lc) {
                if (!lc) return false;
                lc->device = (device == "GPU") ? solver::GPU : solver::CPU;
                lc->method = (method == "INDIRECT") ? solver::INDIRECT : solver::DIRECT;
                return true;
            };

            if (configure(dynamic_cast<loadcase::LinearBuckling*>(base))) return;
            if (configure(dynamic_cast<loadcase::LinearStaticTopo*>(base))) return;
            if (configure(dynamic_cast<loadcase::LinearStatic*>(base))) return;

            throw std::runtime_error("SOLVER not supported for loadcase type " + parser.active_loadcase_type());
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands

