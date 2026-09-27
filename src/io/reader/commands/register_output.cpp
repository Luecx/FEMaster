/**
 * @file register_output.cpp
 * @brief Registers Abaqus/CalculiX-compatible field-output requests.
 *
 * FEMaster intentionally does not distinguish between ASCII and binary output
 * request cards. Abaqus-style *OUTPUT, FIELD with *NODE OUTPUT /
 * *ELEMENT OUTPUT and CalculiX-style *NODE FILE / *EL FILE all modify the same
 * per-step OutputRequestHandler.
 *
 * NODE PRINT and EL PRINT remain outside this result-file interface.
 *
 * @see io::writer::OutputRequestHandler
 * @see io::writer::OutputField
 */

#include "register_functions.h"

#include "../parser.h"
#include "../../dsl/registry.h"
#include "../../writer/output_field.h"

#include "../../../core/logging.h"

#include <array>
#include <string>

namespace fem::io::reader::commands {

namespace {

void request_tokens(Parser& parser, const std::array<std::string, 32>& tokens) {
    auto* loadcase = parser.active_loadcase();
    logging::error(loadcase != nullptr,
        "Output requests require an active analysis step");

    for (const std::string& token : tokens) {
        if (token.empty()) continue;

        const auto field = io::writer::output_field_from_request(token);
        logging::error(field.has_value(),
            "Unsupported output field request: ", token);

        loadcase->output.request(*field);
    }
}

void register_field_list(io::dsl::Registry& registry,
                         Parser& parser,
                         const std::string& command_name,
                         const std::string& description) {
    registry.command(command_name, [&](io::dsl::Command& command) {
        command.allow_if(io::dsl::Condition::parent_is({"LOADCASE", "STEP"}));
        command.doc(description);

        command.variant(io::dsl::Variant::make()
            .segment(io::dsl::Segment::make()
                .range(io::dsl::LineRange{}.min(1))
                .pattern(io::dsl::Pattern::make()
                    .fixed<std::string, 32>().name("FIELD")
                        .desc("Requested Abaqus/CalculiX/FEMaster output variables")
                        .on_missing(std::string{}).on_empty(std::string{})
                )
                .bind([&parser](const std::array<std::string, 32>& fields) {
                    request_tokens(parser, fields);
                })
            )
        );
    });
}

} // anonymous namespace

/**
 * Registers the common result-output grammar.
 *
 * *OUTPUT, FIELD starts an explicit request set in Abaqus syntax. CalculiX
 * cards do not require that container, therefore the first NODE FILE/EL FILE
 * request also starts the explicit set through OutputRequestHandler::request().
 */
void register_output(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("OUTPUT", [&](io::dsl::Command& command) {
        command.allow_if(io::dsl::Condition::parent_is({"LOADCASE", "STEP"}));
        command.doc("Begin an explicit field-output selection for the active step.");

        command.keyword(
            io::dsl::KeywordSpec::make()
                .flag("FIELD").doc("Select field output")
                .key("VARIABLE").optional().allowed({"PRESELECT"})
                    .doc("Restore the step-defined default output selection")
        );

        command.on_enter([&parser](const io::dsl::Keys& keys) {
            auto* loadcase = parser.active_loadcase();
            logging::error(loadcase != nullptr,
                "OUTPUT must follow the analysis procedure that creates the active step");
            logging::error(keys.has("FIELD"),
                "Only *OUTPUT, FIELD is supported");

            if (keys.has("VARIABLE") && keys.raw("VARIABLE") == "PRESELECT") {
                loadcase->output.use_defaults();
            } else {
                loadcase->output.begin_explicit_requests();
            }
        });

        command.variant(io::dsl::Variant::make());
    });

    // Abaqus field-output syntax.
    register_field_list(registry, parser, "NODEOUTPUT",
        "Request nodal result fields for the active step.");
    register_field_list(registry, parser, "ELEMENTOUTPUT",
        "Request element-derived result fields for the active step.");

    // CalculiX result-file syntax. FILE and OUTPUT intentionally have identical
    // semantics in FEMaster; the selected ResultWriters decide the file format.
    register_field_list(registry, parser, "NODEFILE",
        "Request nodal result fields using CalculiX NODE FILE syntax.");
    register_field_list(registry, parser, "ELFILE",
        "Request element-derived result fields using CalculiX EL FILE syntax.");
    register_field_list(registry, parser, "ELOUTPUT",
        "Request element-derived result fields using CalculiX EL OUTPUT syntax.");
}

} // namespace fem::io::reader::commands
