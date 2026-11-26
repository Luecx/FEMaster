/**
 * @file load_c.cpp
 * @brief Implements concentrated nodal force and moment assembly.
 *
 * A concentrated load is already defined in the discrete nodal space and
 * therefore requires no finite-element integration. The implementation replaces
 * omitted `NaN` components by zero, applies the caller's multiplier and
 * transforms local force and moment vectors into the global basis at every
 * target node.
 *
 * With local-to-global basis `A(x_i)` and amplitude `a(t)`, the generalized
 * nodal contribution is
 *
 *     [F_i]       [A(x_i)   0   ] [F_local]
 *     [M_i] += a  [  0    A(x_i)] [M_local].
 *
 * Without an assigned orientation, `A = I`.
 *
 * @see CLoad
 * @see Condition
 * @see cos::CoordinateSystem
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#include "load_c.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>
#include <utility>

namespace fem::bc {

/**
 * Assembles the concentrated generalized load on every target node.
 *
 * The stored six-component vector is split into translational force and nodal
 * moment triplets. `NaN` entries represent omitted components and are replaced
 * by exact zeros before arithmetic is performed. The optional amplitude is
 * replaced by a uniform scale provided by the calling analysis.
 *
 * If an orientation is present, its basis is evaluated independently at every
 * nodal position because curvilinear coordinate systems may vary spatially. A
 * local vector `v_local` is transformed by
 *
 *     v_global = A(x_i) v_local
 *
 * before it is accumulated in the corresponding global generalized DOFs.
 *
 * A dead load evaluates the basis at the reference position. A follower load
 * evaluates it at the current position and contributes the analytic tangent
 *
 *     dF_i/dx_i = d(A(x) F_local)/dx,   dM_i/dx_i = d(A(x) M_local)/dx.
 *
 * Without an orientation the global direction is constant, so a follower
 * load has a vanishing tangent. A follower load without a requested load
 * tangent is reported as an error in either case.
 *
 * @param model_data Reference and current nodal geometry.
 * @param load Generalized nodal load field modified in place.
 * @param load_tangent Optional load-tangent triplets receiving df/du.
 * @param system_dof_ids Global equation numbering of the load tangent.
 * @param scale Multiplier of the nominal load components.
 */
void CLoad::apply(
    model::ModelData&        model_data,
    model::Field&            load,
    TripletList*             load_tangent,
    constraint::Equations&,
    const model::SystemDofs& system_dof_ids,
    TripletList&,
    Precision                scale
) {
    // Validate the definition and the configuration required by the flag
    logging::error(region_ != nullptr,
        "CLOAD: target node region is not initialized");
    logging::error(model_data.positions != nullptr && model_data.positions_reference != nullptr,
        "CLOAD: positions fields are not initialized");
    logging::error(!follower_ || load_tangent != nullptr,
        "CLOAD: follower load requires a load tangent and is only supported in nonlinear analyses");

    // Omitted generalized components contribute no nodal forces or moments
    Vec6 values = Vec6::Zero();
    for (Index i = 0; i < 6; ++i) {
        if (std::isfinite(values_[i])) values[i] = scale * values_[i];
    }

    if (values.isZero()) return;

    // Dead loads use the reference basis, follower loads the current one
    const auto& positions = follower_ ? *model_data.positions : *model_data.positions_reference;

    const Vec3 force  = values.head<3>();
    const Vec3 moment = values.tail<3>();

    for (Index node_id : *region_) {
        // Without an orientation, the global components are constant
        if (!orientation_) {
            for (Index i = 0; i < 6; ++i) {
                load(node_id, i) += values[i];
            }
            continue;
        }

        // Transform the force and moment blocks into global coordinates at the
        // position of the node in the selected configuration
        const Vec3 position    = positions.row_vec3(node_id);
        const Vec3 local_point = orientation_->to_local(position);
        const auto axes        = orientation_->get_axes(local_point);

        const Vec3 global_force  = axes * force;
        const Vec3 global_moment = axes * moment;

        for (Index i = 0; i < 3; ++i) {
            load(node_id, i    ) += global_force [i];
            load(node_id, i + 3) += global_moment[i];
        }

        if (!follower_) continue;

        // The rotating basis couples both blocks to the nodal translation
        const Mat3 force_gradient  = orientation_->get_axes_gradient(position, force);
        const Mat3 moment_gradient = orientation_->get_axes_gradient(position, moment);

        add_load_tangent_block(*load_tangent, system_dof_ids, node_id, 0, node_id, force_gradient);
        add_load_tangent_block(*load_tangent, system_dof_ids, node_id, 3, node_id, moment_gradient);
    }
}


/**
 * Builds a compact diagnostic representation of the concentrated load.
 *
 * The representation reports the semantic node target, all six nominal
 * generalized components and the optional orientation and amplitude. Values are
 * printed before temporal scaling and before local-to-global transformation.
 *
 * @return Human-readable concentrated-load description.
 */
std::string CLoad::str() const {
    std::ostringstream os;

    // Report the stored semantic definition without performing assembly
    os << "CLOAD: target=NSET "
       << (region_ ? region_->name : std::string("?")) << " ("
       << (region_ ? region_->size() : 0) << "), values=["
       << values_[0] << ", "
       << values_[1] << ", "
       << values_[2] << ", "
       << values_[3] << ", "
       << values_[4] << ", "
       << values_[5] << "]";

    if (orientation_) os << ", orientation=" << orientation_->name;
    if (amplitude_  ) os << ", amplitude="   << amplitude_  ->name;

    return os.str();
}

} // namespace fem::bc
