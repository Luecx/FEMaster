// register_loadcase_begin.inl — registers *LOADCASE (begin block)

#include <memory>
#include <stdexcept>
#include <string>

#include "../parser.h"

#include "../../loadcase/linear_buckling.h"
#include "../../loadcase/linear_eigenfreq.h"
#include "../../loadcase/linear_static.h"
#include "../../loadcase/linear_static_topo.h"

namespace fem::input_decks::commands {

inline void register_loadcase_begin(fem::dsl::Registry& registry, Parser& parser) {
    registry.command("LOADCASE", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::parent_is("ROOT"));
        command.doc("Begin a load case definition block.");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("TYPE").required().allowed({
                    "LINEARSTATIC", "LINEARBUCKLING", "LINEARSTATICTOPO", "EIGENFREQ"})
                .key("NAME").optional()
        );

        command.on_enter([&parser](const fem::dsl::Keys& keys) {
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
            } else if (type == "LINEARBUCKLING") {
                lc = std::make_unique<loadcase::LinearBuckling>(id, &wrt, &mdl, 10);
            } else if (type == "LINEARSTATICTOPO") {
                lc = std::make_unique<loadcase::LinearStaticTopo>(id, &wrt, &mdl);
            } else if (type == "EIGENFREQ") {
                lc = std::make_unique<loadcase::LinearEigenfrequency>(id, &wrt, &mdl, 10);
            } else {
                throw std::runtime_error("Unsupported loadcase type: " + type);
            }

            parser.set_active_loadcase(std::move(lc), type);
        });

        command.variant(fem::dsl::Variant::make());
    });
}

} // namespace fem::input_decks::commands

