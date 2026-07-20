/**
 * @file load_c.cpp
 * @brief Implements concentrated nodal force and moment assembly.
 *
 * The implementation converts optional `NaN` component markers into zero,
 * evaluates temporal amplitude scaling, optionally rotates local force and
 * moment triplets at each node, and accumulates the resulting six generalized
 * components in the supplied boundary-condition field.
 *
 * @see load_c.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_c.h"

#include "../core/logging.h"
#include "../model/model_data.h"

#include <cmath>
#include <sstream>
#include <utility>

namespace fem {
namespace bc {
namespace {

/**
 * Converts a partially specified vector into an assembly-safe vector.
 *
 * Input decks use NaN to mark omitted components. Those entries become zero,
 * while the returned flag records whether any component was prescribed.
 *
 * @param vec Vector containing finite values and optional NaN markers.
 * @return The sanitized vector and an active-component flag.
 */
std::pair<Vec3, bool> sanitize_vector(Vec3 vec) {
    bool active = false;

    // Inspect every component independently because forces and moments may be
    // prescribed sparsely in the input deck.
    for (int i = 0; i < 3; ++i) {
        if (std::isnan(vec[i])) {
            // Unspecified components contribute exactly zero during assembly.
            vec[i] = 0.0;
        } else {
            // At least one finite component makes the corresponding triplet
            // relevant for the current node.
            active = true;
        }
    }

    return {vec, active};
}
} // namespace

/**
 * Assembles the concentrated force and moment carried by this load.
 *
 * Sparse components are sanitized, optionally amplitude-scaled and accumulated
 * at every target node. An assigned coordinate system rotates force and moment
 * vectors from its local basis into global coordinates before assembly.
 *
 * @param model_data Model fields and topology required by the load.
 * @param bc Generalized nodal field receiving the contribution.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether amplitude scaling is disabled.
 */
void CLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude) {
    // Positions are required only when a local orientation must be evaluated,
    // but validating them unconditionally keeps the load contract explicit.
    logging::error(model_data.positions != nullptr, "positions field not set in model data");
    logging::error(region_ != nullptr, "CLoad: target node region not set");
    const auto& node_positions = *model_data.positions;

    // Evaluate the shared amplitude once because the same scalar multiplier
    // applies to every node and every active force or moment component.
    const Precision scale = amplitude_ && !ignore_amplitude ? amplitude_->evaluate(time) : Precision(1);

    // The nominal generalized load is repeated for every node in the region.
    for (const ID node_id : *region_) {
        // The nodal position is needed to evaluate position-dependent axes of a
        // local coordinate system.
        const Vec3 position = node_positions.row_vec3(static_cast<Index>(node_id));

        // Split the six generalized components into translational force and
        // rotational moment vectors and replace omitted entries by zero.
        auto [force_local , force_active]  = sanitize_vector(values_.head<3>());
        auto [moment_local, moment_active] = sanitize_vector(values_.tail<3>());

        // Scale only after sanitization so omitted components remain exact zeros
        // even for unusual amplitude values.
        force_local *= scale;
        moment_local *= scale;

        if (!orientation_) {
            // Without a local orientation, the stored components already refer
            // to the global basis and can be accumulated directly.
            if (force_active) {
                for (int i = 0; i < 3; ++i) {
                    bc(node_id, i) += force_local[i];
                }
            }

            if (moment_active) {
                for (int i = 0; i < 3; ++i) {
                    // Rotational components occupy columns three through five
                    // of the generalized nodal field.
                    bc(node_id, static_cast<Dim>(i + 3)) += moment_local[i];
                }
            }

            // The oriented assembly path is unnecessary for this node.
            continue;
        }

        // Map the global nodal position into the coordinate system's local
        // parameterization and evaluate its local basis at that point.
        const Vec3 local_point = orientation_->to_local(position);
        const auto axes        = orientation_->get_axes(local_point);

        if (force_active) {
            // The basis matrix maps local vector components into global
            // Cartesian components before they are added to the translational
            // columns of the load field.
            const Vec3 global_force = axes * force_local;
            for (int i = 0; i < 3; ++i) {
                bc(node_id, i) += global_force[i];
            }
        }

        if (moment_active) {
            // Moments transform with the same orthonormal basis as force
            // vectors, but are accumulated in the rotational columns.
            const Vec3 global_moment = axes * moment_local;
            for (int i = 0; i < 3; ++i) {
                bc(node_id, static_cast<Dim>(i + 3)) += global_moment[i];
            }
        }
    }
}

/**
 * Builds the diagnostic representation of the concentrated load.
 *
 * The result identifies the target node region and reports the nominal force
 * and moment components used to define the load.
 *
 * @return Human-readable load description.
 */
std::string CLoad::str() const {
    std::ostringstream os;

    // Include the target node-set identity, its current cardinality and all six
    // nominal generalized components, including `NaN` placeholders.
    os << "CLOAD: target=NSET "
       << (region_ ? region_->name : std::string("?")) << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0) << "), values=["
       << values_[0] << ", "
       << values_[1] << ", "
       << values_[2] << ", "
       << values_[3] << ", "
       << values_[4] << ", "
       << values_[5] << "]";

    // Append optional modifiers only when they are assigned.
    if (orientation_) {
        os << ", orientation=" << orientation_->name;
    }

    if (amplitude_) {
        os << ", amplitude=" << amplitude_->name;
    }

    return os.str();
}
} // namespace bc
} // namespace fem
