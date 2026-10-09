/**
 * @file load_v.cpp
 * @brief Implements consistent integration of distributed body-force densities.
 *
 * `VLoad` converts a prescribed force density to equivalent nodal forces through
 * the element weak-form contribution
 *
 *     f_e = integral_Omega_e N^T b dOmega.
 *
 * The structural element owns interpolation, quadrature and the physical volume
 * Jacobian. This implementation supplies the vector field `b(x)`, including
 * optional amplitude scaling and local-to-global transformation.
 *
 * The stored vector already has units of force per volume. Consequently the
 * element integration is called with density scaling disabled; multiplying by
 * material density here would incorrectly turn the load into a mass-dependent
 * inertia term.
 *
 * @see VLoad
 * @see Condition
 * @see model::StructuralElement
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#include "load_v.h"

#include "../../core/logging.h"
#include "../../model/element/element_structural.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>

namespace fem::bc {

/**
 * Integrates the prescribed body-force density over all selected elements.
 *
 * Without an amplitude, the effective force density is interpolated between
 * values_start_ and values_ using step_progress. An explicit amplitude instead
 * scales the target at time; ignore_amplitude uses the nominal target. NaN
 * target components are omitted and undefined start components are zero.
 *
 * Each structural element assembles the consistent contribution
 *
 *     f_e = integral_Omega N^T(x) * b(x) dOmega.
 *
 * The integrator evaluates the shape functions, quadrature and physical
 * volume Jacobian. This method provides the global force-density vector
 * b(x), optionally transforming a local vector by the spatially varying
 * coordinate-system basis A(x). Density scaling is disabled because b already
 * represents force per unit volume, rather than acceleration.
 *
 * @param model_data Compiled structural elements and nodal geometry.
 * @param rhs Global nodal force field modified in place.
 * @param equations Constraint equations left unchanged.
 * @param system_dof_ids Global DOF numbering left unchanged.
 * @param matrix Sparse matrix triplets left unchanged.
 * @param time Physical time used for amplitude evaluation.
 * @param ignore_amplitude Assemble the nominal target force density when true.
 * @param step_progress Normalized progress from initial to target load.
 */
void VLoad::apply(
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
        "VLOAD: positions field is not initialized");
    logging::error(region_ != nullptr,
        "VLOAD: target element region is not initialized");

    // Without an amplitude, linearly interpolate between the initial and
    // target force densities. An assigned amplitude scales only the target.
    Precision start_scale = Precision(1) - step_progress;
    Precision end_scale   = step_progress;

    if (ignore_amplitude) {
        start_scale = Precision(0);
        end_scale   = Precision(1);
    } else if (amplitude_) {
        start_scale = Precision(0);
        end_scale   = amplitude_->evaluate(time);
    }

    // Construct the effective force-density vector [bx, by, bz]. Omitted
    // target components contribute zero; undefined initial components are zero.
    Vec3 body_force = Vec3::Zero();

    for (Dim i = 0; i < 3; ++i) {
        if (std::isnan(values_[i])) continue;

        const Precision start = std::isfinite(values_start_[i])
            ? values_start_[i] : Precision(0);

        body_force[i] = start_scale * start + end_scale * values_[i];
    }

    if (body_force.isZero()) return;

    // Convert the distributed body force into equivalent nodal forces:
    //
    //     f_i += sum_q N_i(x_q) * b(x_q) * det(J(x_q)) * w_q
    //
    // The structural integrator evaluates the shape functions N_i and the
    // geometric quadrature measure det(J) * w_q. Density scaling is disabled
    // because body_force already has units of force per unit volume.
    //
    // For an assigned orientation, the callback supplies
    //
    //     b_global(x_q) = A(x_q) * b_local,
    //
    // evaluating the local-to-global axes separately at every quadrature point.
    for (ID element_id : *region_) {
        const auto& element = model_data.elements[element_id];
        if (!element) continue;

        auto* structural = element->as<model::StructuralElement>();
        if (!structural) continue;

        structural->integrate_vector_field(
            rhs,
            false,
            [&](const Vec3& position) -> Vec3 {
                if (!orientation_) return body_force;

                const Vec3 local_point = orientation_->to_local(position);
                const auto axes        = orientation_->get_axes(local_point);

                return axes * body_force;
            }
        );
    }
}

/**
 * Builds a compact diagnostic representation of the body-force density.
 *
 * The output reports the semantic element target, nominal vector and optional
 * orientation and amplitude. The values are shown before coordinate
 * transformation and temporal scaling.
 *
 * @return Human-readable volume-load description.
 */
std::string VLoad::str() const {
    std::ostringstream os;

    // Report the stored force density without scaling or transformation.
    os << "VLOAD: target=ELSET "
       << (region_ ? region_->name : std::string("?")) << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0) << "), values=["
       << values_[0] << ", " << values_[1] << ", " << values_[2] << "]";

    if (orientation_) os << ", orientation=" << orientation_->name;
    if (amplitude_  ) os << ", amplitude="   << amplitude_  ->name;

    return os.str();
}

} // namespace fem::bc
