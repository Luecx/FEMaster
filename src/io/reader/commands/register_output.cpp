/**
 * @file register_output.cpp
 * @brief Registers Abaqus/CalculiX-compatible field-output requests.
 *
 * FEMaster intentionally does not distinguish between ASCII and binary output
 * request cards. NODE OUTPUT / ELEMENT OUTPUT and CalculiX-style NODE FILE /
 * EL FILE all modify the same per-step OutputRequestHandler.
 *
 * Abaqus *OUTPUT, FIELD is accepted as a compatibility container but has no
 * semantic effect. NODE PRINT and EL PRINT are accepted and ignored because
 * FEMaster does not implement a separate text-print result path.
 *
 * @see io::writer::OutputRequestHandler
 * @see io::writer::OutputField
 */

#include "register_functions.h"

#include "../parser.h"
#include "../../dsl/registry.h"
#include "../../writer/output_field.h"

#include "../../../core/logging.h"
#include "../../../loadcase/linear_harmonic.h"

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

        // Harmonic response stores real and imaginary solution components as
        // separate primary fields. Standard Abaqus/CalculiX U/S/E requests
        // therefore expand to both component outputs instead of introducing a
        // synthetic generic displacement/stress/strain source.
        if (loadcase->as<loadcase::LinearHarmonic>() != nullptr) {
            switch (*field) {
                case io::writer::OutputField::DISPLACEMENT:
                    loadcase->output.request(io::writer::OutputField::DISPLACEMENT_REAL);
                    loadcase->output.request(io::writer::OutputField::DISPLACEMENT_IMAG);
                    continue;
                case io::writer::OutputField::STRESS:
                    loadcase->output.request(io::writer::OutputField::STRESS_REAL);
                    loadcase->output.request(io::writer::OutputField::STRESS_IMAG);
                    continue;
                case io::writer::OutputField::STRAIN:
                    loadcase->output.request(io::writer::OutputField::STRAIN_REAL);
                    loadcase->output.request(io::writer::OutputField::STRAIN_IMAG);
                    continue;
                default:
                    break;
            }
        }

        loadcase->output.request(*field);
    }
}

void register_field_list(io::dsl::Registry& registry,
                         Parser& parser,
                         const std::string& command_name,
                         const std::string& description) {
    registry.command(command_name, [&](io::dsl::Command& command) {
        command.allow_if(io::dsl::Condition::parent_is({"LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        command.doc(description);

        command.data(
            io::dsl::Pattern::make()
                .fixed<std::string, 32>("FIELD", "Requested Abaqus/CalculiX/FEMaster output variables")
                    .defaults(std::string{}),
            [&parser](const std::array<std::string, 32>& fields) {
                request_tokens(parser, fields);
            }
        );
    });
}

} // anonymous namespace

/**
 * Registers the common result-output grammar.
 *
 * NODE/ELEMENT result cards are self-contained requests in FEMaster. The first
 * such card replaces the owning step's default request set. Abaqus' surrounding
 * *OUTPUT, FIELD card is therefore syntactic compatibility only and is ignored.
 */
void register_output(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("OUTPUT", [](io::dsl::Command& command) {
        command.allow_if(io::dsl::Condition::parent_is({"LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
        command.doc("Accept the Abaqus OUTPUT container without changing result requests.");

        command.keyword(
            io::dsl::KeywordSpec::make()
                .flag("FIELD").doc("Accepted Abaqus field-output selector")
                .key("VARIABLE").optional().doc("Accepted Abaqus output option")
        );

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

    // Abaqus/CalculiX PRINT output is a separate textual reporting mechanism.
    // FEMaster has no equivalent text-print path, so accept these cards only to
    // keep compatible decks readable and deliberately ignore their variables.
    const auto register_ignored_print = [&](const std::string& name) {
        registry.command(name, [](io::dsl::Command& command) {
            command.allow_if(io::dsl::Condition::parent_is({"LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS"}));
            command.doc("Accept and ignore a textual PRINT output request.");

            command.data(
                io::dsl::Pattern::make()
                    .fixed<std::string, 32>("FIELD")
                        .defaults(std::string{}),
                [](const std::array<std::string, 32>&) {},
                io::dsl::LineRange{}.min(0)
            );
        });
    };

    register_ignored_print("NODEPRINT");
    register_ignored_print("ELPRINT");
}

} // namespace fem::io::reader::commands
