/**
 * @file load_p.cpp
 * @brief Implements consistent pressure integration over finite-element surfaces.
 *
 * The scalar pressure is converted to the geometric traction
 *
 *     t(x) = -p n(x),
 *
 * where `n(x)` is the unit normal computed from the nodal geometry supplied by
 * `ModelData::positions`. Each selected surface then evaluates the consistent
 * weak-form contribution
 *
 *     f_e = -integral_Gamma_e N^T p n dGamma.
 *
 * Surface interpolation, quadrature and the physical surface Jacobian remain
 * encapsulated by the surface implementation. `PLoad` supplies only the scalar
 * magnitude, amplitude scaling and normal traction field.
 *
 * @see PLoad
 * @see Condition
 * @see model::SurfaceInterface
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#include "load_p.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>

namespace fem::bc {

/**
 * Integrates the pressure traction over all selected physical surfaces.
 *
 * Positive pressure acts against the geometric surface normal:
 *
 *     t(x) = -p * n(x).
 *
 * Without an amplitude, pressure is interpolated between its initial and
 * target values using step_progress. An explicit amplitude scales the target
 * pressure at physical time; ignore_amplitude applies the nominal target.
 * An undefined initial pressure is treated as zero.
 *
 * The surface integrator evaluates the consistent nodal force contribution
 *
 *     f_e = integral_Gamma N^T(x) * t(x) dGamma
 *
 * by quadrature. It handles the surface shape functions and area Jacobian;
 * this method supplies the traction vector at each integration point. The
 * unit normal is calculated from model_data.positions, using the same nodal
 * geometry that is passed to the surface integrator.
 *
 * @param model_data Nodal geometry and compiled surface storage.
 * @param rhs Global nodal force field modified in place.
 * @param equations Constraint equations left unchanged.
 * @param system_dof_ids Global DOF numbering left unchanged.
 * @param matrix Sparse matrix triplets left unchanged.
 * @param time Physical time used for amplitude evaluation.
 * @param ignore_amplitude Assemble the nominal target pressure when true.
 * @param step_progress Normalized progress from initial to target pressure.
 */
void PLoad::apply(
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
        "PLOAD: positions field is not initialized");
    logging::error(region_ != nullptr,
        "PLOAD: target surface region is not initialized");
    logging::error(std::isfinite(pressure_),
        "PLOAD: pressure must be finite");

    // Determine the start and target pressure weights. Without an amplitude,
    // p(lambda) = (1 - lambda) * p_start + lambda * p_end.
    // An explicit amplitude replaces this interpolation with p = A(time) * p_end.
    Precision start_scale = Precision(1) - step_progress;
    Precision end_scale   = step_progress;

    if (ignore_amplitude) {
        start_scale = Precision(0);
        end_scale   = Precision(1);
    } else if (amplitude_) {
        start_scale = Precision(0);
        end_scale   = amplitude_->evaluate(time);
    }

    // An undefined initial pressure represents a zero starting load.
    const Precision start = std::isfinite(pressure_start_)
        ? pressure_start_ : Precision(0);
    const Precision pressure = start_scale * start + end_scale * pressure_;

    if (pressure == Precision(0)) return;

    // Convert the scalar pressure into a consistent nodal force for each
    // selected surface:
    //
    //     f_i += -sum_q N_i(x_q) * p * n(x_q) * J_s(x_q) * w_q
    //
    // The surface integrator supplies the shape functions N_i, quadrature
    // weights w_q and physical area Jacobian J_s. The callback returns only
    // the pressure traction -p * n(x_q), using the normal of the provided
    // nodal geometry. Recovering local coordinates is necessary to evaluate
    // the geometric normal at the integration point.
    const auto& positions = *model_data.positions;

    for (ID surface_id : *region_) {
        const auto& surface = model_data.surfaces[surface_id];
        if (!surface) continue;

        surface->integrate_vector_field(
            positions,
            rhs,
            [&](const Vec3& position) -> Vec3 {
                const Vec2 local  = surface->global_to_local(position, positions);
                const Vec3 normal = surface->normal(positions, local);

                return -pressure * normal;
            }
        );
    }
}

/**
 * Builds a compact diagnostic representation of the pressure definition.
 *
 * The output reports the semantic surface target, nominal scalar pressure and
 * optional amplitude. It does not evaluate surface normals or temporal scaling.
 *
 * @return Human-readable pressure-load description.
 */
std::string PLoad::str() const {
    std::ostringstream os;

    // Report the stored pressure without temporal scaling or integration.
    os << "PLOAD: target=SFSET "
       << (region_ ? region_->name : std::string("?")) << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0) << "), p=" << pressure_;

    if (amplitude_) os << ", amplitude=" << amplitude_->name;

    return os.str();
}

} // namespace fem::bc
