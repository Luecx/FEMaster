/**
 * @file load_c.cpp
 * @brief Implements concentrated nodal force and moment assembly.
 *
 * A concentrated load is already defined in the discrete nodal space and
 * therefore requires no finite-element integration. The implementation replaces
 * omitted `NaN` components by zero, evaluates the optional scalar amplitude and
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
 * evaluated once because it scales the complete load uniformly.
 *
 * If an orientation is present, its basis is evaluated independently at every
 * nodal position because curvilinear coordinate systems may vary spatially. A
 * local vector `v_local` is transformed by
 *
 *     v_global = A(x_i) v_local
 *
 * before it is accumulated in the corresponding global generalized DOFs.
 *
 * @param model_data Global nodal geometry used for coordinate transformation.
 * @param rhs Generalized nodal RHS field modified in place.
 * @param equations Constraint-equation output left unchanged by this structural condition.
 * @param system_dof_ids Active equation numbering reserved for condition matrix terms.
 * @param matrix Sparse condition matrix output left unchanged by this structural condition.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Apply the nominal load with unit amplitude when true.
 */
void CLoad::apply(
    model::ModelData&      model_data,
    model::Field&          rhs,
    constraint::Equations&,
    const SystemDofIds&,
    TripletList&,
    Precision              time,
    bool                   ignore_amplitude,
    Precision              step_progress
) {
    logging::error(model_data.positions != nullptr,
        "CLOAD: positions field is not initialized");
    logging::error(region_ != nullptr,
        "CLOAD: target node region is not initialized");

    // Determine the contributions of the start and target values.
    // An explicit amplitude overrides the linear step interpolation.
    Precision start_scale = Precision(1) - step_progress;
    Precision end_scale   = step_progress;

    if (ignore_amplitude) {
        start_scale = Precision(0);
        end_scale   = Precision(1);
    } else if (amplitude_) {
        start_scale = Precision(0);
        end_scale   = amplitude_->evaluate(time);
    }

    // Assemble the local generalized load vector [Fx, Fy, Fz, Mx, My, Mz].
    // NaN target components are omitted; undefined start values are zero.
    Vec6 values = Vec6::Zero();

    // fill in the 6 components if they are not NaN, otherwise they remain zero
    for (Dim i = 0; i < 6; ++i) {
        if (std::isnan(values_[i])) continue;

        const Precision start = std::isfinite(values_start_[i])
            ? values_start_[i] : Precision(0);

        values[i] = start_scale * start + end_scale * values_[i];
    }

    // if no components are active, return early
    if (values.isZero()) return;

    // Add the load to every target node. If an orientation is assigned,
    // transform the force and moment blocks into global coordinates.
    for (ID node_id : *region_) {
        Vec6 global = values;

        if (orientation_) {
            const Vec3 position    = model_data.positions->row_vec3(node_id);
            const Vec3 local_point = orientation_->to_local(position);
            const auto axes        = orientation_->get_axes(local_point);

            global.head<3>() = axes * values.head<3>();
            global.tail<3>() = axes * values.tail<3>();
        }

        for (Dim i = 0; i < 6; ++i) {
            rhs(node_id, i) += global[i];
        }
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
       << (region_ ? static_cast<int>(region_->size()) : 0) << "), values=["
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
