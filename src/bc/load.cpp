/**
 * @file load.cpp
 * @brief Implements the concrete load boundary conditions.
 *
 * Each load type converts high-level user input into equivalent nodal loads
 * or boundary-condition entries by traversing the associated FEM regions.
 *
 * @see src/bc/load.h
 * @see src/bc/load_collector.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load.h"

#include <cmath>
#include <utility>
#include "../model/model_data.h"
#include "../model/element/element.h"
#include "../model/element/element_structural.h"

namespace fem {
namespace bc {

namespace {

std::pair<Vec3, bool> sanitize_vector(Vec3 vec) {
    bool active = false;
    for (int i = 0; i < 3; ++i) {
        if (std::isnan(vec[i])) {
            vec[i] = 0.0;
        } else {
            active = true;
        }
    }
    return {vec, active};
}

} // namespace

/**
 * @copydoc CLoad::apply
 */
void CLoad::apply(model::ModelData& model_data, NodeData& bc, Precision time) {
    auto& node_positions = model_data.get(model::POSITION);
    const Precision scale = amplitude ? amplitude->evaluate(time) : 1.0;

    for (auto& node_id : *region) {
        const Vec3 position = node_positions.row(node_id).head<3>();

        auto [force_local, force_active] = sanitize_vector(values.head<3>());
        auto [moment_local, moment_active] = sanitize_vector(values.tail<3>());

        force_local *= scale;
        moment_local *= scale;

        if (!orientation) {
            if (force_active) {
                for (int i = 0; i < 3; ++i) {
                    bc(node_id, i) += force_local[i];
                }
            }
            if (moment_active) {
                for (int i = 0; i < 3; ++i) {
                    bc(node_id, static_cast<Dim>(i + 3)) += moment_local[i];
                }
            }
            continue;
        }

        const Vec3 local_point = orientation->to_local(position);
        const auto axes = orientation->get_axes(local_point);

        if (force_active) {
            const Vec3 global_force = axes * force_local;
            for (int i = 0; i < 3; ++i) {
                bc(node_id, i) += global_force[i];
            }
        }

        if (moment_active) {
            const Vec3 global_moment = axes * moment_local;
            for (int i = 0; i < 3; ++i) {
                bc(node_id, static_cast<Dim>(i + 3)) += global_moment[i];
            }
        }
    }
}

/**
 * @copydoc DLoad::apply
 */
void DLoad::apply(model::ModelData& model_data, NodeData& bc, Precision time) {
    auto& node_positions = model_data.get(model::POSITION);
    auto [local_values, has_values] = sanitize_vector(values);
    if (!has_values) {
        return;
    }

    const Precision scale = amplitude ? amplitude->evaluate(time) : 1.0;
    local_values *= scale;

    for (auto& surf_id : *region) {
        auto surface = model_data.surfaces[surf_id];
        if (!surface) {
            continue;
        }

        if (!orientation) {
            surface->apply_dload(node_positions, bc, local_values);
            continue;
        }

        const auto contributions = surface->shape_function_integral(node_positions);
        int local_idx = 0;
        for (auto node_it = surface->begin(); node_it != surface->end(); ++node_it, ++local_idx) {
            const ID node_id = *node_it;
            const Vec3 position = node_positions.row(node_id).head<3>();
            const Vec3 local_point = orientation->to_local(position);
            const auto axes = orientation->get_axes(local_point);
            const Vec3 global_values = axes * local_values;

            bc(node_id, 0) += contributions(local_idx) * global_values[0];
            bc(node_id, 1) += contributions(local_idx) * global_values[1];
            bc(node_id, 2) += contributions(local_idx) * global_values[2];
        }
    }
}

/**
 * @copydoc PLoad::apply
 */
void PLoad::apply(model::ModelData& model_data, NodeData& bc, Precision time) {
    auto& node_positions = model_data.get(model::POSITION);
    const Precision scale = amplitude ? amplitude->evaluate(time) : 1.0;
    const Precision scaled_pressure = pressure * scale;
    for (auto& surf_id : *region) {
        model_data.surfaces[surf_id]->apply_pload(node_positions, bc, scaled_pressure);
    }
}

/**
 * @copydoc VLoad::apply
 */
void VLoad::apply(model::ModelData& model_data, NodeData& bc, Precision time) {
    auto& node_positions = model_data.get(model::POSITION);
    auto [local_values, has_values] = sanitize_vector(values);
    if (!has_values) {
        return;
    }

    const Precision scale = amplitude ? amplitude->evaluate(time) : 1.0;
    local_values *= scale;

    for (auto& el_id : *region) {
        auto& el_ptr = model_data.elements[el_id];
        if (!el_ptr) {
            continue;
        }

        auto structural = el_ptr->as<model::StructuralElement>();
        if (!structural) {
            continue;
        }

        Vec3 load_to_apply = local_values;
        if (orientation) {
            Vec3 centroid = Vec3::Zero();
            int node_count = 0;
            for (ID node_id : *structural) {
                centroid += node_positions.row(node_id).head<3>();
                ++node_count;
            }
            if (node_count > 0) {
                centroid /= static_cast<Precision>(node_count);
            }
            const Vec3 local_point = orientation->to_local(centroid);
            const auto axes = orientation->get_axes(local_point);
            load_to_apply = axes * local_values;
        }

        structural->apply_vload(bc, load_to_apply);
    }
}

/**
 * @copydoc TLoad::apply
 */
void TLoad::apply(model::ModelData& model_data, NodeData& bc, Precision time) {
    (void)time;
    for (auto& element_ptr : model_data.elements) {
        if (auto structural = element_ptr->as<model::StructuralElement>()) {
            structural->apply_tload(bc, *temp_field, ref_temp);
        }
    }
}

} // namespace bc
} // namespace fem
