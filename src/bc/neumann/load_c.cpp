/**
 * @file load_c.cpp
 * @brief Implements concentrated nodal force and moment assembly.
 *
 * The implementation converts optional `NaN` component markers into zero,
 * evaluates temporal amplitude scaling, optionally rotates local force and
 * moment triplets at each node, and accumulates the resulting six generalized
 * components in the supplied structural right-hand-side field.
 *
 * @see load_c.h
 * @see StructuralNeumann
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_c.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>
#include <utility>

namespace fem::bc {
namespace {

/**
 * Converts a partially specified vector into an assembly-safe vector.
 *
 * Input definitions use `NaN` to mark omitted components. Those entries become
 * exact zeros so they cannot contaminate subsequent scaling or coordinate
 * transformations, while the returned flag records whether any component was
 * prescribed at all.
 *
 * @param vec Vector containing finite values and optional NaN markers.
 * @return Sanitized vector and flag indicating at least one active component.
 */
std::pair<Vec3, bool> sanitize_vector(Vec3 vec) {
    bool active = false;

    // Inspect every component independently because force and moment triplets may
    // be prescribed sparsely in the input definition.
    for (int i = 0; i < 3; ++i) {
        if (std::isnan(vec[i])) {
            // Unspecified components contribute exactly zero during assembly.
            vec[i] = Precision(0);
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
 * vectors from its local basis into global coordinates before assembly. The
 * target must be a NODE field with the six generalized structural components.
 *
 * @param model_data Model fields and topology required by the load.
 * @param rhs Generalized nodal right-hand-side field receiving the contribution.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether amplitude scaling is disabled.
 */
void CLoad::apply(model::ModelData& model_data,
                  model::Field&     rhs,
                  Precision         time,
                  bool              ignore_amplitude) {
    // Validate all shared model and target data before evaluating any node.
    logging::error(model_data.positions != nullptr,
        "positions field not set in model data");
    logging::error(region_ != nullptr,
        "CLoad: target node region not set");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components >= 6,
        "CLoad: target field must be a NODE field with at least six components");

    const auto& node_positions = *model_data.positions;

    // Evaluate the common amplitude once because the same scalar multiplier
    // applies to every target node and every active generalized component.
    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);

    // Apply the nominal generalized load independently to every node in the
    // selected region.
    for (const ID node_id : *region_) {
        // The physical node position is required only when a position-dependent
        // local coordinate system must be evaluated.
        const Vec3 position = node_positions.row_vec3(static_cast<Index>(node_id));

        // Split the six generalized components into force and moment triplets and
        // replace omitted entries by exact zeros before any arithmetic.
        auto [force_local,  force_active]  = sanitize_vector(values_.head<3>());
        auto [moment_local, moment_active] = sanitize_vector(values_.tail<3>());

        // Scale only the sanitized vectors so NaN sentinels never enter the
        // numerical right-hand-side field.
        force_local  *= scale;
        moment_local *= scale;

        if (!orientation_) {
            // Without a local orientation the stored components already refer to
            // the global basis and can be accumulated directly.
            if (force_active) {
                for (Dim i = 0; i < 3; ++i) rhs(node_id, i) += force_local[i];
            }

            if (moment_active) {
                // Rotational components occupy columns three through five of the
                // generalized structural nodal field.
                for (Dim i = 0; i < 3; ++i) rhs(node_id, i + 3) += moment_local[i];
            }
            continue;
        }

        // Evaluate the potentially position-dependent local basis at the current
        // node before transforming the prescribed vector components.
        const Vec3 local_point = orientation_->to_local(position);
        const auto axes        = orientation_->get_axes(local_point);

        if (force_active) {
            // Map local force components into global Cartesian coordinates and
            // accumulate them in the translational structural components.
            const Vec3 global_force = axes * force_local;
            for (Dim i = 0; i < 3; ++i) rhs(node_id, i) += global_force[i];
        }

        if (moment_active) {
            // Moments transform with the same local basis but are written to the
            // rotational generalized components.
            const Vec3 global_moment = axes * moment_local;
            for (Dim i = 0; i < 3; ++i) rhs(node_id, i + 3) += global_moment[i];
        }
    }
}

/**
 * Builds the diagnostic representation of the concentrated load.
 *
 * The result identifies the target node region and reports the nominal force and
 * moment components before amplitude scaling or coordinate transformation.
 *
 * @return Human-readable load description.
 */
std::string CLoad::str() const {
    std::ostringstream os;

    // Include the target node-set identity, its current cardinality and all six
    // nominal generalized components, including NaN placeholders.
    os << "CLOAD: target=NSET "
       << (region_ ? region_->name : std::string("?")) << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0) << "), values=["
       << values_[0] << ", " << values_[1] << ", " << values_[2] << ", "
       << values_[3] << ", " << values_[4] << ", " << values_[5] << "]";

    // Append optional modifiers only when they are assigned.
    if (orientation_) os << ", orientation=" << orientation_->name;
    if (amplitude_) os << ", amplitude=" << amplitude_->name;

    return os.str();
}

} // namespace fem::bc
