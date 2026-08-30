/**
 * @file load_c.cpp
 * @brief Implements concentrated nodal force and moment assembly.
 *
 * The implementation converts optional `NaN` component markers into zero,
 * evaluates temporal amplitude scaling, optionally rotates local force and
 * moment triplets at each node, and accumulates the resulting six generalized
 * components in the supplied right-hand-side field.
 *
 * `CLoad` never modifies the system operator. The optional DOF map and LHS
 * triplet list supplied by the common `Neumann::apply()` contract are therefore
 * intentionally ignored.
 *
 * @see load_c.h
 * @author Finn Eggers
 * @date 30.08.2026
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
 * Input decks use NaN to mark omitted components. Those entries become zero,
 * while the returned flag records whether any component was prescribed.
 *
 * @param vec Vector containing finite values and optional NaN markers.
 * @return Sanitized vector and an active-component flag.
 */
std::pair<Vec3, bool> sanitize_vector(Vec3 vec) {
    bool active = false;

    for (int i = 0; i < 3; ++i) {
        if (std::isnan(vec[i])) {
            vec[i] = Precision(0);
        } else {
            active = true;
        }
    }

    return {vec, active};
}

} // namespace

/**
 * Assembles the concentrated force and moment carried by this condition.
 *
 * Sparse components are sanitized, optionally amplitude-scaled and accumulated
 * at every target node. An assigned coordinate system rotates force and moment
 * vectors from its local basis into global coordinates before assembly.
 *
 * @param model_data Model fields and topology required by the condition.
 * @param rhs Six-component structural nodal RHS receiving the contribution.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether amplitude scaling is disabled.
 * @param system_dof_ids Unused optional system DOF map.
 * @param lhs Unused optional system-matrix triplet list.
 */
void CLoad::apply(model::ModelData&       model_data,
                  model::Field&           rhs,
                  Precision               time,
                  bool                    ignore_amplitude,
                  const SystemDofIds*      system_dof_ids,
                  TripletList*             lhs) {
    (void) system_dof_ids;
    (void) lhs;

    // Validate the nodal target and position field used by optional orientations.
    logging::error(model_data.positions != nullptr,
        "positions field not set in model data");
    logging::error(region_ != nullptr,
        "CLoad: target node region not set");
    const auto& node_positions = *model_data.positions;

    // Evaluate the shared amplitude once for the complete condition.
    const Precision scale = amplitude_ && !ignore_amplitude ? amplitude_->evaluate(time) : Precision(1);

    for (const ID node_id : *region_) {
        const Vec3 position = node_positions.row_vec3(static_cast<Index>(node_id));

        // Split generalized components and replace omitted entries by zero.
        auto [force_local , force_active ] = sanitize_vector(values_.head<3>());
        auto [moment_local, moment_active] = sanitize_vector(values_.tail<3>());
        force_local  *= scale;
        moment_local *= scale;

        if (!orientation_) {
            if (force_active) {
                for (Dim i = 0; i < 3; ++i) {
                    rhs(node_id, i) += force_local[i];
                }
            }
            if (moment_active) {
                for (Dim i = 0; i < 3; ++i) {
                    rhs(node_id, i + 3) += moment_local[i];
                }
            }
            continue;
        }

        // Evaluate the local basis at the current nodal position.
        const Vec3 local_point = orientation_->to_local(position);
        const auto axes        = orientation_->get_axes(local_point);

        if (force_active) {
            const Vec3 global_force = axes * force_local;
            for (Dim i = 0; i < 3; ++i) {
                rhs(node_id, i) += global_force[i];
            }
        }
        if (moment_active) {
            const Vec3 global_moment = axes * moment_local;
            for (Dim i = 0; i < 3; ++i) {
                rhs(node_id, i + 3) += global_moment[i];
            }
        }
    }
}

/**
 * Builds the diagnostic representation of the concentrated load.
 *
 * @return Human-readable load description.
 */
std::string CLoad::str() const {
    std::ostringstream os;
    os << "CLOAD: target=NSET "
       << (region_ ? region_->name : std::string("?")) << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0) << "), values=["
       << values_[0] << ", " << values_[1] << ", " << values_[2] << ", "
       << values_[3] << ", " << values_[4] << ", " << values_[5] << "]";

    if (orientation_)
        os << ", orientation=" << orientation_->name;
    if (amplitude_)
        os << ", amplitude=" << amplitude_->name;

    return os.str();
}

} // namespace fem::bc
