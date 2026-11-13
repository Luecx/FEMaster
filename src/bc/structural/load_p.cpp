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
 * The calling analysis supplies a multiplier for the nominal pressure.
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
 */
void PLoad::apply(
    model::ModelData&      model_data,
    model::Field&          rhs,
    constraint::Equations&,
    const SystemDofIds&    system_dof_ids,
    TripletList&           matrix,
    Precision              scale
) {
    logging::error(model_data.positions != nullptr,
        "PLOAD: positions field is not initialized");
    logging::error(region_ != nullptr,
        "PLOAD: target surface region is not initialized");
    logging::error(std::isfinite(pressure_),
        "PLOAD: pressure must be finite");

    const Precision pressure = scale * pressure_;

    if (pressure == Precision(0)) return;

    // Consistent normal traction on the reference or current geometry.
    // The surface routine computes the oriented area vector and its exact
    // displacement derivative in a single quadrature pass.
    const auto& positions = follower_
        ? *model_data.positions : *model_data.positions_reference;
    const SystemDofIds empty_ids;
    const auto& tangent_ids = follower_ ? system_dof_ids : empty_ids;
    for (ID surface_id : *region_) {
        const auto& surface = model_data.surfaces[surface_id];
        if (!surface) continue;
        surface->integrate_pressure(positions, rhs, pressure, tangent_ids, matrix);
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
