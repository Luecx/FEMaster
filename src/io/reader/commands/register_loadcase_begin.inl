/**
 * @file register_loadcase_begin.inl
 * @brief Registers creation and completion of FEMaster load-case scopes.
 *
 * The root-level `LOADCASE` command maps the requested analysis type to one
 * concrete `loadcase::LoadCase` implementation and transfers its ownership to
 * `Parser::begin_loadcase()`. Consecutive child commands then configure that
 * active object through the common parser lifecycle.
 *
 * When the scope closes, `Parser::end_loadcase()` removes the active definition,
 * executes the selected analysis and releases its storage. This file therefore
 * owns the common construction boundary, while formulation-specific settings
 * remain in their dedicated registration files.
 *
 * @see loadcase::LoadCase
 * @see Parser::begin_loadcase
 * @see Parser::end_loadcase
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include <memory>
#include <string>

#include "../parser.h"

#include "../../../core/logging.h"
#include "../../../loadcase/linear_buckling.h"
#include "../../../loadcase/linear_eigenfreq.h"
#include "../../../loadcase/linear_harmonic.h"
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
                    "LINEARHARMONIC", "NONLINEARSTATIC"})
                .key("NAME").optional()
        );

        command.on_enter([&parser](const fem::io::dsl::Keys& keys) {
            logging::error(parser.active_loadcase() == nullptr,
                "Nested *LOADCASE blocks are not supported");

            const std::string type = keys.raw("TYPE");
            if (type == "LINEARSTATIC") {
                parser.begin_loadcase(std::make_unique<loadcase::LinearStatic>());
            } else if (type == "NONLINEARSTATIC") {
                parser.begin_loadcase(std::make_unique<loadcase::NonlinearStatic>());
            } else if (type == "LINEARBUCKLING") {
                parser.begin_loadcase(std::make_unique<loadcase::LinearBuckling>());
            } else if (type == "LINEARSTATICTOPO") {
                parser.begin_loadcase(std::make_unique<loadcase::LinearStaticTopo>());
            } else if (type == "EIGENFREQ") {
                parser.begin_loadcase(std::make_unique<loadcase::LinearEigenfrequency>());
            } else if (type == "LINEARTRANSIENT") {
                parser.begin_loadcase(std::make_unique<loadcase::Transient>());
            } else if (type == "LINEARHARMONIC") {
                parser.begin_loadcase(std::make_unique<loadcase::LinearHarmonic>());
            } else {
                logging::error(false,
                    "Unsupported loadcase type: ", type);
            }
        });

        // Execute and release the completed load case when its scope closes
        command.on_exit([&parser](const fem::io::dsl::Keys&) {
            if (!parser.active_loadcase()) {
                // Ignore a load case already consumed by an earlier closing command
                return;
            }
            parser.end_loadcase();
        });

        command.variant(fem::io::dsl::Variant::make());
    });
}

} // namespace fem::io::reader::commands
