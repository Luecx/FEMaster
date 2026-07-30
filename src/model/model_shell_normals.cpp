/**
 * @file model_shell_normals.cpp
 * @brief Applies explicit element-nodal fields to shell reference normals.
 */

#include "model.h"

#include "../core/logging.h"

#include <cmath>
#include <vector>

namespace fem::model {

void Model::apply_shell_element_normal_field(const Field::Ptr& field,
                                             Precision equalize_angle_degrees) {
    logging::error(field != nullptr, "Model: shell normal field is null");
    logging::error(field->domain == FieldDomain::ELEMENT_NODAL,
                   "Model: shell normal field '", field->name,
                   "' must use ELEMENT_NODAL domain");
    logging::error(field->components == 3,
                   "Model: shell normal field '", field->name,
                   "' must have exactly three components");
    logging::error(_data->shell_element_nodal_normals != nullptr,
                   "Model: automatic shell normals have not been built");
    logging::error(field->rows == _data->shell_element_nodal_normals->rows,
                   "Model: shell normal field '", field->name,
                   "' has an incompatible row count");

    constexpr Precision normal_tolerance = Precision(1e-12);

    const Precision pi             = std::acos(Precision(-1));
    const Precision equalize_angle = equalize_angle_degrees * pi / Precision(180);
    const Precision cos_equalize   = std::cos(equalize_angle);

    const Field::Ptr automatic = _data->shell_element_nodal_normals;

    std::vector<std::vector<Index>> node_rows(static_cast<std::size_t>(_data->max_nodes));
    std::vector<Vec3> row_normals(static_cast<std::size_t>(field->rows), Vec3::Zero());
    std::vector<bool> explicit_rows(static_cast<std::size_t>(field->rows), false);

    // Read valid prescribed vectors and use the automatically generated normal
    // as fallback for every unset, non-finite or nearly zero field row.
    for (const ElementPtr& element : _data->elements) {
        if (!element) {
            continue;
        }

        const auto structural = element->as<StructuralElement>();
        if (!structural || !structural->is_shell()) {
            continue;
        }

        const Index offset  = static_cast<Index>(structural->elem_nodal_offset);
        const Index n_nodes = static_cast<Index>(structural->n_nodes());

        for (Index local_node = 0; local_node < n_nodes; ++local_node) {
            const ID node_id = structural->nodes()[local_node];
            const Index row  = offset + local_node;

            Vec3 normal = field->row_vec3(row);
            const bool prescribed = normal.allFinite() && normal.norm() > normal_tolerance;

            if (prescribed) {
                normal.normalize();
                explicit_rows[static_cast<std::size_t>(row)] = true;
            } else {
                normal = automatic->row_vec3(row);
                logging::error(normal.allFinite() && normal.norm() > normal_tolerance,
                               "Model: invalid automatic shell normal in element ",
                               structural->elem_id, " at local node ", local_node);
                normal.normalize();
            }

            row_normals[static_cast<std::size_t>(row)] = normal;
            node_rows[static_cast<std::size_t>(node_id)].push_back(row);
        }
    }

    // Use the same angular clustering as the automatic normal construction.
    // Explicit rows remain fixed. Automatic rows follow the explicit average
    // inside their cluster, or the ordinary cluster average when no override
    // is present.
    for (const auto& rows : node_rows) {
        std::vector<Vec3> cluster_normals;
        std::vector<Precision> cluster_weights;
        std::vector<std::vector<Index>> cluster_rows;

        for (Index row : rows) {
            const Vec3 normal = row_normals[static_cast<std::size_t>(row)];
            bool added = false;

            for (Index cluster = 0; cluster < static_cast<Index>(cluster_normals.size()); ++cluster) {
                if (normal.dot(cluster_normals[static_cast<std::size_t>(cluster)]) < cos_equalize) {
                    continue;
                }

                Vec3& cluster_normal = cluster_normals[static_cast<std::size_t>(cluster)];
                Precision& weight   = cluster_weights[static_cast<std::size_t>(cluster)];

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
            Vec3 explicit_target = Vec3::Zero();
            bool has_explicit = false;

            for (Index row : cluster_rows[static_cast<std::size_t>(cluster)]) {
                if (!explicit_rows[static_cast<std::size_t>(row)]) {
                    continue;
                }
                explicit_target += row_normals[static_cast<std::size_t>(row)];
                has_explicit = true;
            }

            Vec3 target = cluster_normals[static_cast<std::size_t>(cluster)];
            if (has_explicit) {
                logging::error(explicit_target.norm() > normal_tolerance,
                               "Model: prescribed shell normals cancel inside one angular cluster");
                target = explicit_target.normalized();
            }

            for (Index row : cluster_rows[static_cast<std::size_t>(cluster)]) {
                if (!explicit_rows[static_cast<std::size_t>(row)]) {
                    row_normals[static_cast<std::size_t>(row)] = target;
                }
            }
        }
    }

    // Store both prescribed and automatically completed values in the selected
    // generic field and expose it as the semantic shell-normal field.
    for (Index row = 0; row < field->rows; ++row) {
        const Vec3& normal = row_normals[static_cast<std::size_t>(row)];
        if (normal.norm() <= normal_tolerance) {
            continue;
        }
        (*field)(row, 0) = normal(0);
        (*field)(row, 1) = normal(1);
        (*field)(row, 2) = normal(2);
    }

    _data->shell_element_nodal_normals = field;
}

} // namespace fem::model
