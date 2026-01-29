// register_field.inl — registers *FIELD

#include <array>
#include <sstream>
#include <limits>
#include <stdexcept>
#include <string>

#include "../../core/logging.h"
#include "../../core/types_num.h"
#include "../../data/field.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"
#include "../../model/model.h"

namespace fem::input_decks::commands {

namespace detail {

constexpr std::size_t kMaxFieldCols = 64;

inline model::FieldDomain parse_field_domain(const std::string& raw) {
    if (raw == "NODE") return model::FieldDomain::NODE;
    if (raw == "ELEMENT") return model::FieldDomain::ELEMENT;
    if (raw == "IP") return model::FieldDomain::IP;
    throw std::runtime_error("FIELD: unknown TYPE '" + raw + "' (expected NODE/ELEMENT/IP)");
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
    if (ss.fail() || !ss.eof()) {
        throw std::runtime_error("FIELD: failed to parse value '" + token + "'");
    }
    return value;
}

} // namespace detail

inline void register_field(fem::dsl::Registry& registry, model::Model& model) {
    registry.command("FIELD", [&](fem::dsl::Command& command) {
        command.allow_if(fem::dsl::Condition::any_of({
            fem::dsl::Condition::parent_is("ROOT"),
            fem::dsl::Condition::parent_is("LOADCASE")
        }));
        command.doc("Create or populate a generic field (NODE/ELEMENT/IP).");

        command.keyword(
            fem::dsl::KeywordSpec::make()
                .key("NAME").required().doc("Field name")
                .key("TYPE").required().allowed({"NODE","ELEMENT","IP"}).doc("Field domain")
                .key("COLS").required().doc("Number of components")
                .key("FILL").optional("ZERO").allowed({"ZERO","NAN"}).doc("Initialization (ZERO default)")
        );

        struct Context {
            model::Field::Ptr field = nullptr;
            Index cols = 0;
        };
        auto ctx = std::make_shared<Context>();

        command.on_enter([&model, ctx](const fem::dsl::Keys& keys) {
            const std::string name = keys.raw("NAME");
            const std::string type = keys.raw("TYPE");
            const std::string fill = keys.has("FILL") ? keys.raw("FILL") : std::string("ZERO");
            const int cols = keys.get<int>("COLS");

            if (cols <= 0) {
                throw std::runtime_error("FIELD: COLS must be > 0");
            }
            if (static_cast<std::size_t>(cols) > detail::kMaxFieldCols) {
                throw std::runtime_error("FIELD: COLS exceeds internal limit of " + std::to_string(detail::kMaxFieldCols));
            }

            const auto domain = detail::parse_field_domain(type);
            const bool fill_nan = (fill == "NAN");

            if (domain == model::FieldDomain::IP && model._data->max_integration_points == 0) {
                throw std::runtime_error("FIELD: TYPE=IP requires max_integration_points to be configured");
            }

            ctx->field = model._data->create_field(name, domain, static_cast<Index>(cols), fill_nan);
            ctx->cols = static_cast<Index>(cols);

            if (fill_nan) {
                ctx->field->fill_nan();
            } else {
                ctx->field->set_zero();
            }
        });

        static const std::string kSkipToken = "__SKIP__";

        command.variant(fem::dsl::Variant::make()
            .segment(fem::dsl::Segment::make()
                .range(fem::dsl::LineRange{}.min(0))
                .pattern(fem::dsl::Pattern::make()
                    .one<fem::ID>().name("ID").desc("Row id")
                    .fixed<std::string, detail::kMaxFieldCols>().name("V").desc("Values")
                        .on_empty(kSkipToken).on_missing(kSkipToken)
                )
                .bind([ctx](fem::ID id, const std::array<std::string, detail::kMaxFieldCols>& values) {
                    if (!ctx->field) {
                        throw std::runtime_error("FIELD: internal error (field not initialized)");
                    }

                    if (id < 0 || static_cast<Index>(id) >= ctx->field->rows) {
                        logging::error(false, "FIELD: id ", id, " out of bounds for field '",
                                       ctx->field->name, "' (rows=", ctx->field->rows, ")");
                        return;
                    }

                    const Index row = static_cast<Index>(id);
                    for (Index c = 0; c < ctx->cols; ++c) {
                        const auto& token = values[static_cast<std::size_t>(c)];
                        if (token == kSkipToken) {
                            continue;
                        }
                        const Precision v = detail::parse_precision_or_throw(token);
                        (*ctx->field)(row, c) = v;
                    }

                    for (Index c = ctx->cols; c < static_cast<Index>(detail::kMaxFieldCols); ++c) {
                        if (values[static_cast<std::size_t>(c)] != kSkipToken) {
                            throw std::runtime_error("FIELD: more values than COLS for field '" + ctx->field->name + "'");
                        }
                    }
                })
            )
        );
    });
}

} // namespace fem::input_decks::commands
