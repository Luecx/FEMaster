/**
 * @file register_initial_condition.cpp
 * @brief Registers model-level initial temperature and velocity fields.
 *
 * The reader assigns either an existing named node field or independent node
 * storage populated directly by this command's data rows. Inline fields carry
 * diagnostic names but are not published in the ModelData field registry.
 * Model compilation supplies node counts and instance-qualified row mappings;
 * analyses consume the assigned model state when initializing their solvers.
 *
 * @see model::ModelData
 * @see register_initial_condition
 * @author Finn Eggers
 * @date 11.10.2026
 */

#include "register_functions.h"
#include "../../dsl/registry.h"

#include "../../../core/logging.h"
#include "../../../data/field.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

#include <array>
#include <limits>
#include <memory>
#include <sstream>
#include <string>

namespace fem::io::reader::commands {

/**
 * Registers initial node temperatures and generalized node velocities.
 *
 * FIELD selects existing named storage and admits no data rows. Without FIELD,
 * command entry creates unregistered, zero-initialized node storage; each data
 * row supplies a compiled node reference followed by component values. Node
 * references may be numeric or instance-qualified, as in the NODE form of FIELD.
 * Temperature has one component; velocity uses [vx,vy,vz,wx,wy,wz] in the
 * model's global generalized-DOF basis.
 *
 * Empty and omitted values preserve the current entry, initially zero. Repeated
 * rows overwrite only supplied components. Explicit NaN and infinity use the
 * same parsing semantics as FIELD. Extra components and invalid node references
 * are rejected before row assignment. ModelData retains shared ownership through
 * its temperature or velocity pointer; no named field is created in inline mode.
 * Both command aliases execute after model compilation in the root scope.
 *
 * @param registry Reader registry receiving both initial-condition aliases.
 * @param model Compiled model receiving the initial state field.
 */
void register_initial_condition(fem::io::dsl::Registry& registry, model::Model& model) {
    const auto register_name = [&](const char* name) {
        registry.command(name, [&](fem::io::dsl::Command& command) {
            command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
            command.doc("Assign a named FIELD or read node reference and initial values directly. "
                        "TEMPERATURE uses one component; VELOCITY uses vx,vy,vz,wx,wy,wz. "
                        "Inline fields start at zero; empty or omitted values preserve entries.");

            command.keyword(
                fem::io::dsl::KeywordSpec::make()
                    .key("TYPE").required().allowed({"TEMPERATURE", "VELOCITY"})
                    .key("FIELD").optional()
            );

            // Keep the current field alive across command entry and data callbacks.
            // Each occurrence replaces this handle; ModelData retains assigned fields.
            auto field = std::make_shared<model::Field::Ptr>();

            command.on_enter([&model, field](const fem::io::dsl::Keys& keys) {
                const std::string type       = keys.raw("TYPE");
                const std::string field_name = keys.raw("FIELD");
                const Index       components = type == "TEMPERATURE" ? 1 : 6;

                // Resolve shared named storage or allocate independent zeroed storage.
                // Diagnostic names do not enter the registry because reg is false.
                if (keys.has("FIELD")) {
                    logging::error(!field_name.empty(),
                        "INITIALCONDITION: FIELD must specify a non-empty field name");
                    *field = model._data->get_field(field_name);
                    logging::error(*field != nullptr,
                        "INITIALCONDITION: field ", field_name, " does not exist");
                    logging::error((*field)->domain == model::FieldDomain::NODE,
                        "INITIALCONDITION: field ", field_name, " must use NODE domain");
                    logging::error((*field)->components == components,
                        "INITIALCONDITION: ", type, " field must have ", components, " components");
                } else {
                    *field = model._data->create_field("INITIAL_" + type,
                        model::FieldDomain::NODE, components, false, false);
                }

                // Assign the same storage that subsequent inline data rows populate.
                if (type == "TEMPERATURE") {
                    model._data->temperature = *field;
                } else {
                    model._data->velocity = *field;
                }
            });

            // Referenced fields are assigned as-is and accept no following data rows.
            command.variant(fem::io::dsl::Variant::make()
                .when(fem::io::dsl::Condition::key_present("FIELD")));

            // Inline rows follow FIELD's node addressing and skipped-value semantics.
            const std::string skip_token = "__SKIP__";
            command.variant(fem::io::dsl::Variant::make()
                .when(fem::io::dsl::Condition::negate(fem::io::dsl::Condition::key_present("FIELD")))
                .segment(fem::io::dsl::Segment::make()
                    .range(fem::io::dsl::LineRange{}.min(1))
                    .pattern(fem::io::dsl::Pattern::make()
                        .one<std::string>().name("TARGET").desc("Node ID or INSTANCE.ID")
                        .fixed<std::string, 6>().name("V").desc("Initial component values")
                            .on_empty(skip_token).on_missing(skip_token)
                    )
                    .bind([&model, field, skip_token](const std::string& target,
                                                     const std::array<std::string, 6>& values) {
                        // Resolve the node and validate the active field layout first.
                        logging::error(*field != nullptr,
                            "INITIALCONDITION: field not initialized");
                        const ID node = model.compiled_node_id(target);
                        logging::error(node >= 0 && static_cast<Index>(node) < (*field)->rows,
                            "INITIALCONDITION: node ", target, " is out of bounds");
                        for (Index component = (*field)->components; component < values.size(); ++component) {
                            logging::error(values[component] == skip_token,
                                "INITIALCONDITION: too many components for ", (*field)->name);
                        }

                        // Parse supplied components locally, preserving FIELD's explicit
                        // NaN/Inf spellings and leaving skipped components unchanged.
                        for (Index component = 0; component < (*field)->components; ++component) {
                            const std::string& token = values[component];
                            if (token == skip_token) continue;

                            Precision value{};
                            if (token == "NAN" || token == "+NAN" || token == "-NAN") {
                                value = std::numeric_limits<Precision>::quiet_NaN();
                            } else if (token == "INF" || token == "+INF") {
                                value = std::numeric_limits<Precision>::infinity();
                            } else if (token == "-INF") {
                                value = -std::numeric_limits<Precision>::infinity();
                            } else {
                                std::istringstream stream(token);
                                stream >> value;
                                logging::error(!stream.fail() && stream.eof(),
                                    "INITIALCONDITION: failed to parse value '", token, "'");
                            }

                            (**field)(static_cast<Index>(node), component) = value;
                        }
                    })
                )
            );
        });
    };

    register_name("INITIALCONDITION");
    register_name("INITIALCONDITIONS");
}

} // namespace fem::io::reader::commands
