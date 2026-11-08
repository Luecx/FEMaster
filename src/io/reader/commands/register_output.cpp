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

void request_tokens(Parser&                           parser,
                    const std::array<std::string, 32>& tokens,
                    model::FieldDomain                 domain) {
    auto* loadcase = parser.active_loadcase();
    logging::error(loadcase != nullptr,
        "Output requests require an active analysis step");

    for (const std::string& token : tokens) {
        if (token.empty()) continue;

        const auto field = io::writer::output_field_from_request(token);
        if (!field.has_value()) {
            logging::warning(false, "Unsupported output field request: ", token, " (ignored)");
            continue;
        }

        // Use the native element domain for element-wide scalar quantities
        model::FieldDomain requested_domain = domain;
        if (domain == model::FieldDomain::ELEMENT_NODAL) {
            switch (*field) {
                case io::writer::OutputField::COMPLIANCE:
                case io::writer::OutputField::VOLUME:
                case io::writer::OutputField::DENS_GRAD:
                case io::writer::OutputField::DENSITY:
                case io::writer::OutputField::ORIENTATION:
                case io::writer::OutputField::ORIENTATION_GRAD:
                    requested_domain = model::FieldDomain::ELEMENT;
                    break;
                default:
                    break;
            }
        }

        // Non-spatial results have no nodal or element association.
        // LAMBDA remains a frame field; spectral summaries are step metadata.
        if (*field == io::writer::OutputField::LAMBDA
            || io::writer::output_field_is_step_field(*field)) {
            requested_domain = model::FieldDomain::UNKNOWN;
        }

        // Harmonic response stores real and imaginary solution components as
        // separate primary fields. Standard Abaqus/CalculiX U/S/E requests
        // therefore expand to both component outputs instead of introducing a
        // synthetic generic displacement/stress/strain source.
        if (loadcase->as<loadcase::LinearHarmonic>() != nullptr) {
            switch (*field) {
                case io::writer::OutputField::DISPLACEMENT:
                    loadcase->output.request(io::writer::OutputField::DISPLACEMENT_REAL, requested_domain);
                    loadcase->output.request(io::writer::OutputField::DISPLACEMENT_IMAG, requested_domain);
                    continue;
                case io::writer::OutputField::STRESS:
                    loadcase->output.request(io::writer::OutputField::STRESS_REAL, requested_domain);
                    loadcase->output.request(io::writer::OutputField::STRESS_IMAG, requested_domain);
                    continue;
                case io::writer::OutputField::STRAIN:
                    loadcase->output.request(io::writer::OutputField::STRAIN_REAL, requested_domain);
                    loadcase->output.request(io::writer::OutputField::STRAIN_IMAG, requested_domain);
                    continue;
                default:
                    break;
            }
        }

        loadcase->output.request(*field, requested_domain);
    }
}

void register_field_list(io::dsl::Registry& registry,
                         Parser& parser,
                         const std::string& command_name,
                         const std::string& description) {
    const bool element_output = command_name == "ELEMENTOUTPUT"
                             || command_name == "ELFILE"
                             || command_name == "ELOUTPUT";
    const model::FieldDomain domain = element_output
        ? model::FieldDomain::ELEMENT_NODAL
        : model::FieldDomain::NODE;

    registry.command(command_name, [&](io::dsl::Command& command) {
        command.allow_if(io::dsl::Condition::parent_is({"LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS", "HEATTRANSFER"}));
        command.doc(description);

        command.variant(io::dsl::Variant::make()
            .segment(io::dsl::Segment::make()
                .range(io::dsl::LineRange{}.min(1))
                .pattern(io::dsl::Pattern::make()
                    .fixed<std::string, 32>().name("FIELD")
                        .desc("Requested Abaqus/CalculiX/FEMaster output variables")
                        .on_missing(std::string{}).on_empty(std::string{})
                )
                .bind([&parser, domain](const std::array<std::string, 32>& fields) {
                    request_tokens(parser, fields, domain);
                })
            )
        );
    });
}

} // anonymous namespace

/**
 * Registers the common result-output grammar.
 *
 * NODE/ELEMENT result cards are self-contained requests in FEMaster. The first
 * such card appends an explicit request to the owning step. Load-case defaults
 * are added only at execution time when the request list is still empty.
 */
void register_output(fem::io::dsl::Registry& registry, Parser& parser) {
    registry.command("OUTPUT", [](io::dsl::Command& command) {
        command.allow_if(io::dsl::Condition::parent_is({"LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS", "HEATTRANSFER"}));
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
            command.allow_if(io::dsl::Condition::parent_is({"LOADCASE", "STATIC", "FREQUENCY", "BUCKLE", "DYNAMIC", "STEADYSTATEDYNAMICS", "HEATTRANSFER"}));
            command.doc("Accept and ignore a textual PRINT output request.");

            command.variant(io::dsl::Variant::make()
                .segment(io::dsl::Segment::make()
                    .range(io::dsl::LineRange{}.min(0))
                    .pattern(io::dsl::Pattern::make()
                        .fixed<std::string, 32>().name("FIELD")
                            .on_missing(std::string{}).on_empty(std::string{})
                    )
                    .bind([](const std::array<std::string, 32>&) {})
                )
            );
        });
    };

    register_ignored_print("NODEPRINT");
    register_ignored_print("ELPRINT");
}

} // namespace fem::io::reader::commands
