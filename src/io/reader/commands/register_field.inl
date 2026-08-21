/**
 * @file register_field.inl
 * @brief Registers generic fields on the compiled finite-element model.
 *
 * `FIELD` creates named rectangular storage on node, element, surface or
 * element-nodal domains after `Model::compile()` has established deterministic
 * row counts. Data rows may use direct compiled identifiers or qualified deck
 * references that are resolved through the common reader helpers.
 *
 * The command validates component counts and row addressing before writing
 * values into `model::Field`. Registration in `ModelData` makes the completed
 * field available to loads, topology analyses and result-processing routines.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include <array>
#include <limits>
#include <memory>
#include <sstream>
#include <string>

#include "../reference.h"
#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../data/field.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../../model/model.h"

namespace fem::io::reader::commands {

namespace detail {

constexpr std::size_t kMaxFieldCols = 64;

inline model::FieldDomain parse_field_domain(const std::string& raw) {
    if (raw == "NODE") return model::FieldDomain::NODE;
    if (raw == "ELEMENT") return model::FieldDomain::ELEMENT;
    if (raw == "ELEMENT_NODAL" || raw == "ELEMENTNODAL") return model::FieldDomain::ELEMENT_NODAL;
    if (raw == "ELEMENT_IP" || raw == "ELEMENTIP" || raw == "IP") return model::FieldDomain::ELEMENT_IP;
    if (raw == "ELEMENT_MP" || raw == "ELEMENTMP" || raw == "MP") return model::FieldDomain::ELEMENT_MP;

    logging::error(false,
        "FIELD: unknown TYPE '", raw,
        "' (expected NODE/ELEMENT/ELEMENT_NODAL/ELEMENT_IP/ELEMENT_MP)");
    return model::FieldDomain::NODE;
}

inline Precision parse_precision_or_throw(const std::string& token) {
    if (token == "NAN" || token == "+NAN" || token == "-NAN") {
        return std::numeric_limits<Precision>::quiet_NaN();
    }
    if (token == "INF" || token == "+INF") {
        return std::numeric_limits<Precision>::infinity();
    }
    if (token == "-INF") {
        return -std::numeric_limits<Precision>::infinity();
    }

    std::istringstream ss(token);
    Precision value{};
    ss >> value;
    logging::error(!ss.fail() && ss.eof(),
        "FIELD: failed to parse value '", token, "'");
    return value;
}

} // namespace detail

inline void register_field(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("FIELD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::any_of({
            fem::io::dsl::Condition::parent_is("ROOT"),
            fem::io::dsl::Condition::parent_is("LOADCASE")
        }));
        command.doc("Create or populate a generic field using semantic instance-local entity references.");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").required().doc("Field name")
                .key("TYPE").required()
                    .allowed({"NODE", "ELEMENT", "ELEMENT_NODAL", "ELEMENTNODAL",
                              "ELEMENT_IP", "ELEMENTIP", "IP", "ELEMENT_MP", "ELEMENTMP", "MP"})
                    .doc("Field domain")
                .key("COLS").required().doc("Number of components")
                .key("FILL").optional("ZERO").allowed({"ZERO", "NAN"}).doc("Initialization")
        );

        struct Context {
            model::Field::Ptr field = nullptr;
            Index cols = 0;
        };
        auto ctx = std::make_shared<Context>();

        command.on_enter([&model, ctx](const fem::io::dsl::Keys& keys) {
            const std::string name = keys.raw("NAME");
            const std::string type = keys.raw("TYPE");
            const std::string fill = keys.raw("FILL");
            const int cols = keys.get<int>("COLS");

            logging::error(cols > 0,
                "FIELD: COLS must be > 0");
            logging::error(static_cast<std::size_t>(cols) <= detail::kMaxFieldCols,
                "FIELD: COLS exceeds internal limit of ", detail::kMaxFieldCols);

            const auto domain = detail::parse_field_domain(type);
            if (domain == model::FieldDomain::ELEMENT_NODAL) {
                logging::error(model._data->element_nodal_offsets != nullptr,
                    "FIELD: TYPE=ELEMENT_NODAL requires compiled element enumeration");
            }
            if (domain == model::FieldDomain::ELEMENT_IP) {
                logging::error(model._data->element_ip_offsets != nullptr,
                    "FIELD: TYPE=ELEMENT_IP requires compiled integration-point enumeration");
            }
            if (domain == model::FieldDomain::ELEMENT_MP) {
                logging::error(model._data->element_mp_offsets != nullptr,
                    "FIELD: TYPE=ELEMENT_MP requires compiled material-point enumeration");
            }

            const bool fill_nan = fill == "NAN";
            ctx->field = model._data->create_field(name, domain, static_cast<Index>(cols), fill_nan);
            ctx->cols  = static_cast<Index>(cols);

            if (fill_nan) ctx->field->fill_nan();
            else ctx->field->set_zero();
        });

        static const std::string kSkipToken = "__SKIP__";

        const auto assign_values = [ctx](Index row,
                                         const std::array<std::string, detail::kMaxFieldCols>& values) {
            logging::error(ctx->field != nullptr,
                "FIELD: internal error (field not initialized)");
            logging::error(row < ctx->field->rows,
                "FIELD: resolved row ", row, " is out of bounds for field ", ctx->field->name);

            for (Index c = 0; c < ctx->cols; ++c) {
                const auto& token = values[static_cast<std::size_t>(c)];
                if (token != kSkipToken) {
                    (*ctx->field)(row, c) = detail::parse_precision_or_throw(token);
                }
            }
            for (Index c = ctx->cols; c < static_cast<Index>(detail::kMaxFieldCols); ++c) {
                logging::error(values[static_cast<std::size_t>(c)] == kSkipToken,
                    "FIELD: more values than COLS for field ", ctx->field->name);
            }
        };

        const auto get_element = [&model](const std::string& target) {
            const ID element_id = io::reader::compiled_element_id(model, target);
            logging::error(element_id >= 0
                        && static_cast<std::size_t>(element_id) < model._data->elements.size()
                        && model._data->elements[static_cast<std::size_t>(element_id)] != nullptr,
                "FIELD: element ", target, " is not configured");
            return model._data->elements[static_cast<std::size_t>(element_id)];
        };

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"NODE", "ELEMENT"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("TARGET").desc("ID or INSTANCE.ID")
                    .fixed<std::string, detail::kMaxFieldCols>().name("V").desc("Values")
                        .on_empty(kSkipToken).on_missing(kSkipToken)
                )
                .bind([&model, ctx, assign_values](const std::string& target,
                                                   const std::array<std::string, detail::kMaxFieldCols>& values) {
                    const ID row = ctx->field->domain == model::FieldDomain::NODE
                        ? io::reader::compiled_node_id(model, target)
                        : io::reader::compiled_element_id(model, target);
                    assign_values(static_cast<Index>(row), values);
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"ELEMENT_NODAL", "ELEMENTNODAL"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("ELEMENT").desc("Element ID or INSTANCE.ID")
                    .one<ID>().name("LOCAL_NODE").desc("Zero-based local node index")
                    .fixed<std::string, detail::kMaxFieldCols>().name("V").desc("Values")
                        .on_empty(kSkipToken).on_missing(kSkipToken)
                )
                .bind([get_element, assign_values](const std::string& target,
                                                   ID local_node,
                                                   const std::array<std::string, detail::kMaxFieldCols>& values) {
                    const auto element = get_element(target);
                    logging::error(local_node >= 0 && local_node < element->n_nodes(),
                        "FIELD: local node ", local_node, " is out of bounds for element ", target);
                    assign_values(static_cast<Index>(element->elem_nodal_offset + local_node), values);
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"ELEMENT_IP", "ELEMENTIP", "IP"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("ELEMENT").desc("Element ID or INSTANCE.ID")
                    .one<ID>().name("LOCAL_IP").desc("Zero-based local integration-point index")
                    .fixed<std::string, detail::kMaxFieldCols>().name("V").desc("Values")
                        .on_empty(kSkipToken).on_missing(kSkipToken)
                )
                .bind([get_element, assign_values](const std::string& target,
                                                   ID local_ip,
                                                   const std::array<std::string, detail::kMaxFieldCols>& values) {
                    const auto element = get_element(target);
                    logging::error(local_ip >= 0 && local_ip < element->num_ip(),
                        "FIELD: local integration point ", local_ip,
                        " is out of bounds for element ", target);
                    assign_values(element->ip_index(static_cast<Index>(local_ip)), values);
                })
            )
        );

        command.variant(fem::io::dsl::Variant::make()
            .when(fem::io::dsl::Condition::key_equals("TYPE", {"ELEMENT_MP", "ELEMENTMP", "MP"}))
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
                    .one<std::string>().name("ELEMENT").desc("Element ID or INSTANCE.ID")
                    .one<ID>().name("LOCAL_IP").desc("Zero-based local integration-point index")
                    .one<ID>().name("LOCAL_MP").desc("Zero-based local material-point index")
                    .fixed<std::string, detail::kMaxFieldCols>().name("V").desc("Values")
                        .on_empty(kSkipToken).on_missing(kSkipToken)
                )
                .bind([get_element, assign_values](const std::string& target,
                                                   ID local_ip,
                                                   ID local_mp,
                                                   const std::array<std::string, detail::kMaxFieldCols>& values) {
                    const auto element = get_element(target);
                    logging::error(local_ip >= 0 && local_ip < element->num_ip(),
                        "FIELD: local integration point ", local_ip,
                        " is out of bounds for element ", target);
                    logging::error(local_mp >= 0 && static_cast<Index>(local_mp) < element->num_mp_per_ip(),
                        "FIELD: local material point ", local_mp,
                        " is out of bounds for element ", target);
                    assign_values(
                        element->mp_index(static_cast<Index>(local_ip), static_cast<Index>(local_mp)),
                        values
                    );
                })
            )
        );
    });
}

} // namespace fem::io::reader::commands
