// register_loadcase_solver.inl — registers SOLVER within *LOADCASE

#include <stdexcept>
#include <string>

#include "../parser.h"

#include "../../../loadcase/linear_buckling.h"
#include "../../../loadcase/linear_static.h"
#include "../../../loadcase/linear_static_topo.h"
#include "../../../loadcase/linear_transient.h"
#include "../../../loadcase/nonlinear_static.h"

namespace fem::io::reader::commands {

inline void register_loadcase_solver(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("SOLVER", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("LOADCASE"));
        command.doc(
            "Configure solver options for the active loadcase.\n"
            "\n"
            "Constraint | Backend   | DIRECT       | INDIRECT\n"
            "NULLSPACE  | CPU MKL   | Yes          | Yes\n"
            "NULLSPACE  | CPU Eigen | Yes          | Yes\n"
            "NULLSPACE  | GPU       | Yes          | Yes\n"
            "NULLSPACE  | GPU cuDSS | Yes          | Yes\n"
            "LAGRANGE   | CPU MKL   | Yes          | No\n"
            "LAGRANGE   | CPU Eigen | Limited      | No\n"
            "LAGRANGE   | GPU       | No           | No\n"
            "LAGRANGE   | GPU cuDSS | Yes          | No"
        );

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("DEVICE").optional("CPU").allowed({"CPU", "GPU"})
                .key("METHOD").optional("DIRECT").allowed({"DIRECT", "INDIRECT"})
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            auto* base = parser.active_loadcase();
            logging::error(base != nullptr, "SOLVER must appear inside *LOADCASE");

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
            if (configure(dynamic_cast<loadcase::NonlinearStatic*>(base))) return;
            if (configure(dynamic_cast<loadcase::Transient*>(base))) return;

            throw std::runtime_error("SOLVER not supported for loadcase type " + parser.active_loadcase_type());
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
