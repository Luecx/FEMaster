// register_field.inl — registers *FIELD and *NORMAL

#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../../core/logging.h"
#include "../../../core/types_num.h"
#include "../../../data/field.h"
#include "../../../model/element/element_structural.h"
#include "../../../model/model.h"
#include "../../dsl/condition.h"
#include "../../dsl/keyword.h"

namespace fem::io::reader::commands {

namespace detail {

constexpr std::size_t kMaxFieldCols = 64;

inline model::FieldDomain parse_field_domain(const std::string& raw) {
    if (raw == "NODE") return model::FieldDomain::NODE;
    if (raw == "ELEMENT") return model::FieldDomain::ELEMENT;
    if (raw == "ELEMENT_NODAL" || raw == "ELEMENTNODAL") return model::FieldDomain::ELEMENT_NODAL;
    if (raw == "ELEMENT_IP" || raw == "ELEMENTIP" || raw == "IP") return model::FieldDomain::ELEMENT_IP;
    if (raw == "ELEMENT_MP" || raw == "ELEMENTMP" || raw == "MP") return model::FieldDomain::ELEMENT_MP;
    throw std::runtime_error("FIELD: unknown TYPE '" + raw + "' (expected NODE/ELEMENT/ELEMENT_NODAL/ELEMENT_IP/ELEMENT_MP)");
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

inline void register_field(fem::io::dsl::Registry& registry, model::Model& model) {
    registry.command("FIELD", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::any_of({
            fem::io::dsl::Condition::parent_is("ROOT"),
            fem::io::dsl::Condition::parent_is("LOADCASE")
        }));
        command.doc("Create or populate a generic field (NODE/ELEMENT/ELEMENT_NODAL/ELEMENT_IP/ELEMENT_MP).");

        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("NAME").required().doc("Field name")
                .key("TYPE").required().allowed({"NODE","ELEMENT","ELEMENT_NODAL","ELEMENT_IP","IP","ELEMENT_MP","MP"}).doc("Field domain")
                .key("COLS").required().doc("Number of components")
                .key("FILL").optional("ZERO").allowed({"ZERO","NAN"}).doc("Initialization (ZERO default)")
        );

        struct Context {
            model::Field::Ptr field = nullptr;
            Index cols = 0;
        };
        auto ctx = std::make_shared<Context>();

        command.on_enter([&model, ctx](const fem::io::dsl::Keys& keys) {
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

            if (domain == model::FieldDomain::ELEMENT_NODAL &&
                (!model._data->element_nodal_offsets ||
                 (*model._data->element_nodal_offsets)(static_cast<Index>(model._data->max_elems)) <= 0)) {
                throw std::runtime_error("FIELD: TYPE=ELEMENT_NODAL requires elements to be configured");
            }
            if (domain == model::FieldDomain::ELEMENT_IP &&
                (!model._data->element_ip_offsets ||
                 (*model._data->element_ip_offsets)(static_cast<Index>(model._data->max_elems)) <= 0)) {
                throw std::runtime_error("FIELD: TYPE=ELEMENT_IP requires integration points to be configured");
            }
            if (domain == model::FieldDomain::ELEMENT_MP &&
                (!model._data->element_mp_offsets ||
                 (*model._data->element_mp_offsets)(static_cast<Index>(model._data->max_elems)) <= 0)) {
                throw std::runtime_error("FIELD: TYPE=ELEMENT_MP requires material points to be configured");
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

        command.variant(fem::io::dsl::Variant::make()
            .segment(fem::io::dsl::Segment::make()
                .range(fem::io::dsl::LineRange{}.min(0))
                .pattern(fem::io::dsl::Pattern::make()
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
                        (*ctx->field)(row, c) = detail::parse_precision_or_throw(token);
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

    registry.command("NORMAL", [&](fem::io::dsl::Command& command) {
        command.allow_if(fem::io::dsl::Condition::parent_is("ROOT"));
        command.doc("Use an ELEMENT_NODAL vector field as shell reference normals.");
        command.keyword(
            fem::io::dsl::KeywordSpec::make()
                .key("FIELD").required().doc("ELEMENT_NODAL field containing three normal components")
        );
        command.variant(fem::io::dsl::Variant::make());

        command.on_enter([&model](const fem::io::dsl::Keys& keys) {
            constexpr Precision normal_tolerance = Precision(1e-12);

            const std::string field_name = keys.raw("FIELD");
            const auto normals = model._data->get_field(field_name);
            const auto automatic_normals = model._data->shell_element_nodal_normals;

            logging::error(normals != nullptr,
                           "NORMAL: field '", field_name, "' does not exist");
            logging::error(normals->domain == model::FieldDomain::ELEMENT_NODAL,
                           "NORMAL: field '", field_name, "' must use ELEMENT_NODAL domain");
            logging::error(normals->components == 3,
                           "NORMAL: field '", field_name, "' must have exactly three components");
            logging::error(automatic_normals != nullptr,
                           "NORMAL: automatic shell normals are not initialized");
            logging::error(normals->rows == automatic_normals->rows,
                           "NORMAL: field '", field_name, "' has an invalid row count");

            const Precision pi             = std::acos(Precision(-1));
            const Precision equalize_angle = Precision(20) * pi / Precision(180);
            const Precision cos_equalize   = std::cos(equalize_angle);

            std::vector<std::vector<Index>> node_rows(static_cast<std::size_t>(model._data->max_nodes));
            std::vector<Vec3> row_normals(static_cast<std::size_t>(normals->rows), Vec3::Zero());
            std::vector<bool> explicit_normal(static_cast<std::size_t>(normals->rows), false);

            for (const model::ElementPtr& element: model._data->elements) {
                if (!element) {
                    continue;
                }

                const auto structural = element->as<model::StructuralElement>();
                if (!structural || !structural->is_shell()) {
                    continue;
                }

                const Index offset  = static_cast<Index>(structural->elem_nodal_offset);
                const Index n_nodes = static_cast<Index>(structural->n_nodes());

                for (Index local_node = 0; local_node < n_nodes; ++local_node) {
                    const ID node_id = structural->nodes()[local_node];
                    const Index row  = offset + local_node;
                    Vec3 normal      = normals->row_vec3(row);

                    const bool is_explicit = normal.allFinite() && normal.norm() > normal_tolerance;
                    if (!is_explicit) {
                        normal = automatic_normals->row_vec3(row);
                    }

                    logging::error(normal.allFinite() && normal.norm() > normal_tolerance,
                                   "NORMAL: invalid shell normal for element ", structural->elem_id,
                                   " at local node ", local_node);

                    normal.normalize();
                    row_normals[static_cast<std::size_t>(row)] = normal;
                    explicit_normal[static_cast<std::size_t>(row)] = is_explicit;
                    node_rows[static_cast<std::size_t>(node_id)].push_back(row);
                }
            }

            for (const auto& rows: node_rows) {
                std::vector<Vec3> cluster_normals;
                std::vector<Precision> cluster_weights;
                std::vector<std::vector<Index>> cluster_rows;

                for (Index row: rows) {
                    const Vec3 normal = row_normals[static_cast<std::size_t>(row)];
                    bool added = false;

                    for (Index cluster = 0; cluster < static_cast<Index>(cluster_normals.size()); ++cluster) {
                        if (normal.dot(cluster_normals[static_cast<std::size_t>(cluster)]) < cos_equalize) {
                            continue;
                        }

                        Vec3& cluster_normal = cluster_normals[static_cast<std::size_t>(cluster)];
                        Precision& weight = cluster_weights[static_cast<std::size_t>(cluster)];
                        cluster_normal = (weight * cluster_normal + normal).normalized();
                        weight += Precision(1);
                        cluster_rows[static_cast<std::size_t>(cluster)].push_back(row);
                        added = true;
                        break;
                    }

                    if (!added) {
                        cluster_normals.push_back(normal);
                        cluster_weights.push_back(Precision(1));
                        cluster_rows.push_back({row});
                    }
                }

                for (Index cluster = 0; cluster < static_cast<Index>(cluster_rows.size()); ++cluster) {
                    Vec3 explicit_sum = Vec3::Zero();
                    Index explicit_count = 0;

                    for (Index row: cluster_rows[static_cast<std::size_t>(cluster)]) {
                        if (explicit_normal[static_cast<std::size_t>(row)]) {
                            explicit_sum += row_normals[static_cast<std::size_t>(row)];
                            ++explicit_count;
                        }
                    }

                    const Vec3 target = explicit_count > 0
                        ? explicit_sum.normalized()
                        : cluster_normals[static_cast<std::size_t>(cluster)];

                    for (Index row: cluster_rows[static_cast<std::size_t>(cluster)]) {
                        const Vec3 normal = explicit_normal[static_cast<std::size_t>(row)]
                            ? row_normals[static_cast<std::size_t>(row)]
                            : target;
                        (*normals)(row, 0) = normal(0);
                        (*normals)(row, 1) = normal(1);
                        (*normals)(row, 2) = normal(2);
                    }
                }
            }

            model._data->shell_element_nodal_normals = normals;
        });
    });
}

} // namespace fem::io::reader::commands
