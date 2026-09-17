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
 * @see Neumann
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
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Apply the nominal load with unit amplitude when true.
 */
void CLoad::apply(model::ModelData& model_data, model::Field& rhs, Precision time, bool ignore_amplitude) {
    // Validate the target and geometry required by the assembly operation
    logging::error(model_data.positions != nullptr,
        "CLOAD: positions field is not initialized");
    logging::error(region_ != nullptr,
        "CLOAD: target node region is not initialized");

    const auto& node_positions = *model_data.positions;

    // Convert the sparse input convention into a numerical vector. The lambda
    // remains local because this NaN convention is needed only by this assembly
    // routine.
    auto sanitize_vector = [](Vec3 vector) {
        bool active = false;
        for (Dim component = 0; component < 3; ++component) {
            if (std::isnan(vector[component])) {
                vector[component] = Precision(0);
            } else {
                active = true;
            }
        }
        return std::pair<Vec3, bool>{vector, active};
    };

    // Split the generalized load into the two physical vector blocks
    auto [force_local,  force_active ] = sanitize_vector(values_.head<3>());
    auto [moment_local, moment_active] = sanitize_vector(values_.tail<3>());

    if (!force_active && !moment_active) {
        return;
    }

    // Apply the common scalar time history once to the nominal vectors
    const Precision scale = amplitude_ && !ignore_amplitude ? amplitude_->evaluate(time) : Precision(1);
    force_local  *= scale;
    moment_local *= scale;

    // Assemble the same nominal concentrated load on every node in the region
    for (ID node_id : *region_) {
        if (!orientation_) {
            // Global components already coincide with the solver basis
            if (force_active) {
                for (Dim component = 0; component < 3; ++component) {
                    rhs(node_id, component) += force_local[component];
                }
            }

            if (moment_active) {
                for (Dim component = 0; component < 3; ++component) {
                    rhs(node_id, static_cast<Dim>(component + 3)) += moment_local[component];
                }
            }
            continue;
        }

        // Evaluate the local basis A(x_i) at the current nodal position
        const Vec3 position    = node_positions.row_vec3(static_cast<Index>(node_id));
        const Vec3 local_point = orientation_->to_local(position);
        const auto axes        = orientation_->get_axes(local_point);

        // Transform and assemble the translational force block
        if (force_active) {
            const Vec3 force_global = axes * force_local;
            for (Dim component = 0; component < 3; ++component) {
                rhs(node_id, component) += force_global[component];
            }
        }

        // Transform and assemble the nodal moment block with the same basis
        if (moment_active) {
            const Vec3 moment_global = axes * moment_local;
            for (Dim component = 0; component < 3; ++component) {
                rhs(node_id, static_cast<Dim>(component + 3)) += moment_global[component];
            }
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

    if (orientation_) {
        os << ", orientation=" << orientation_->name;
    }

    if (amplitude_) {
        os << ", amplitude=" << amplitude_->name;
    }

    return os.str();
}

} // namespace fem::bc
