// register_loadcase_begin.inl — registers *LOADCASE (begin block)

#include <memory>
#include <stdexcept>
#include <string>

#include "../parser.h"

#include "../../../loadcase/linear_buckling.h"
#include "../../../loadcase/linear_eigenfreq.h"
#include "../../../loadcase/linear_static.h"
#include "../../../loadcase/linear_static_topo.h"
#include "../../../loadcase/linear_transient.h"
#include "../../../loadcase/nonlinear_static.h"

namespace fem::io::reader::commands {

inline void register_loadcase_begin(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("LOADCASE", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Begin a load case definition block.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("TYPE").required().allowed({
                    "LINEARSTATIC", "LINEARBUCKLING", "LINEARSTATICTOPO", "EIGENFREQ", "LINEARTRANSIENT",
                    "NONLINEARSTATIC"})
                .key("NAME").optional()
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            if (parser.active_loadcase()) {
                throw std::runtime_error("Nested *LOADCASE blocks are not supported");
            }

            auto& mdl = parser.model();
            auto& wrt = parser.writer();
            const std::string type = keys.raw("TYPE");
            const int id = parser.next_loadcase_id();

            std::unique_ptr<loadcase::LoadCase> lc;
            if (type == "LINEARSTATIC") {
                lc = std::make_unique<loadcase::LinearStatic>(id, &wrt, &mdl);
            } else if (type == "NONLINEARSTATIC") {
                lc = std::make_unique<loadcase::NonlinearStatic>(id, &wrt, &mdl);
            } else if (type == "LINEARBUCKLING") {
                lc = std::make_unique<loadcase::LinearBuckling>(id, &wrt, &mdl, 10);
            } else if (type == "LINEARSTATICTOPO") {
                lc = std::make_unique<loadcase::LinearStaticTopo>(id, &wrt, &mdl);
            } else if (type == "EIGENFREQ") {
                lc = std::make_unique<loadcase::LinearEigenfrequency>(id, &wrt, &mdl, 10);
            } else if (type == "LINEARTRANSIENT") {
                lc = std::make_unique<loadcase::Transient>(id, &wrt, &mdl);
            } else {
                throw std::runtime_error("Unsupported loadcase type: " + type);
            }

            parser.set_active_loadcase(std::move(lc), type);
        });

        // When the LOADCASE scope exits, run and clear it.
        command.on_exit([&parser](const fem::io::dsl::Keys&) {
            auto* lc = parser.active_loadcase();
            if (!lc) {
                // Nothing to do (either already cleared, or mis-scoped END).
                return;
            }
            try {
                lc->run();
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("LOADCASE execution failed: ") + e.what());
            }
            parser.clear_active_loadcase();
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
