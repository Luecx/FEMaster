/**
 * @file load_d.cpp
 * @brief Implements consistent integration of distributed surface tractions.
 *
 * For every selected surface, the continuous traction field is converted to
 * equivalent nodal forces through the standard weak-form contribution
 *
 *     f_e = integral_Gamma_e N^T t dGamma.
 *
 * The surface implementation owns interpolation, quadrature and the physical
 * area Jacobian. `DLoad` provides the traction integrand, including optional
 * amplitude scaling and local-to-global coordinate transformation. Curvilinear
 * coordinate systems are evaluated independently at every integration point.
 *
 * @see DLoad
 * @see Condition
 * @see model::SurfaceInterface
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#include "load_d.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>

namespace fem::bc {

/**
 * Integrates the effective traction over all selected physical surfaces.
 *
 * The stored target traction is interpolated with the start value using
 * step_progress, unless a named amplitude scales the target at time.
 * ignore_amplitude bypasses both operations and assembles the nominal target.
 * NaN target components are omitted; undefined start components are zero.
 *
 * Each surface integrates N^T t over its physical area. For a local traction,
 * the coordinate-system axes are evaluated at each quadrature point before
 * the integrator accumulates equivalent global nodal forces.
 *
 * @param model_data Global nodal positions and compiled surface geometry.
 * @param rhs Global nodal force field modified in place.
 * @param equations Constraint equations left unchanged.
 * @param system_dof_ids Global DOF numbering left unchanged.
 * @param matrix Sparse matrix triplets left unchanged.
 * @param time Physical time used for amplitude evaluation.
 * @param ignore_amplitude Assemble the nominal target traction when true.
 * @param step_progress Normalized progress between start and target tractions.
 */
void DLoad::apply(
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
        "DLOAD: positions field is not initialized");
    logging::error(region_ != nullptr,
        "DLOAD: target surface region is not initialized");

    // Determine the interpolation weights for the prescribed traction.
    //
    // Without an amplitude, the load changes linearly from its initial
    // value t_start to the target value t_end:
    //
    //     t(lambda) = (1 - lambda) * t_start + lambda * t_end
    //
    // An explicit amplitude instead scales the target traction directly:
    //
    //     t(time) = A(time) * t_end
    //
    // ignore_amplitude bypasses both and uses the nominal target traction.
    Precision start_scale = Precision(1) - step_progress;
    Precision end_scale   = step_progress;

    if (ignore_amplitude) {
        start_scale = Precision(0);
        end_scale   = Precision(1);
    } else if (amplitude_) {
        start_scale = Precision(0);
        end_scale   = amplitude_->evaluate(time);
    }

    // Construct the effective traction vector in the prescribed coordinate
    // system. NaN marks an omitted target component, while an undefined
    // initial component is interpreted as zero.
    Vec3 traction = Vec3::Zero();

    for (Dim i = 0; i < 3; ++i) {
        if (std::isnan(values_[i])) continue;

        const Precision start = std::isfinite(values_start_[i])
            ? values_start_[i] : Precision(0);

        traction[i] = start_scale * start + end_scale * values_[i];
    }

    if (traction.isZero()) return;

    // Convert the distributed surface traction into equivalent nodal forces.
    // For each surface element, the consistent load vector is
    //
    //     f_e = integral_Gamma N^T(x) * t(x) dGamma
    //
    // where N contains the surface shape functions and t(x) is the traction
    // vector expressed in global coordinates.
    //
    // The surface integration routine evaluates this integral numerically:
    //
    //     f_i += sum_q N_i(x_q) * t(x_q) * J_s(x_q) * w_q
    //
    // with N_i the shape function of node i, J_s the surface Jacobian and
    // w_q the quadrature weight. The resulting nodal forces are accumulated
    // directly in rhs.
    //
    // Only the traction field t(x) is supplied here. If orientation_ is set,
    // the local traction is transformed into global coordinates using
    //
    //     t_global(x) = A(x) * t_local
    //
    // The transformation is evaluated at every quadrature point because the
    // local coordinate axes may vary across the surface.
    const auto& positions = *model_data.positions;

    for (ID surface_id : *region_) {
        const auto& surface = model_data.surfaces[surface_id];
        if (!surface) continue;

        surface->integrate_vector_field(
            positions,
            rhs,
            [&](const Vec3& position) -> Vec3 {
                if (!orientation_) return traction;

                const Vec3 local_point = orientation_->to_local(position);
                const auto axes        = orientation_->get_axes(local_point);

                return axes * traction;
            }
        );
    }
}

/**
 * Builds a compact diagnostic representation of the distributed traction.
 *
 * The output reports the semantic surface target, the nominal traction vector
 * before amplitude scaling and coordinate transformation, and the optional
 * modifier names.
 *
 * @return Human-readable traction-load description.
 */
std::string DLoad::str() const {
    std::ostringstream os;

    // Report the stored traction before amplitude scaling or transformation.
    os << "DLOAD: target=SFSET "
       << (region_ ? region_->name : std::string("?")) << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0) << "), values=["
       << values_[0] << ", " << values_[1] << ", " << values_[2] << "]";

    if (orientation_) os << ", orientation=" << orientation_->name;
    if (amplitude_  ) os << ", amplitude="   << amplitude_  ->name;

    return os.str();
}

} // namespace fem::bc
