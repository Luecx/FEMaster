/**
 * @file load_p.cpp
 * @brief Implements consistent pressure integration over finite-element surfaces.
 *
 * The scalar pressure is converted to the geometric traction
 *
 *     t(x) = -p n(x),
 *
 * where `n(x)` is the unit normal of the loaded configuration. Each selected
 * surface then evaluates the consistent weak-form contribution
 *
 *     f_e = -integral_Gamma_e N^T p n dGamma.
 *
 * Dead pressure uses the reference positions. Follower pressure uses the
 * current positions and contributes the analytic load tangent df/du, which the
 * nonlinear analysis subtracts from the structural tangent.
 *
 * Surface interpolation, quadrature, the physical surface Jacobian and its
 * linearization remain encapsulated by the surface implementation. `PLoad`
 * supplies only the scalar magnitude, the configuration and the global
 * scattering.
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
 * The calling analysis supplies a multiplier for the nominal pressure. Each
 * surface evaluates the consistent nodal force
 *
 *     f_e = -p integral_Gamma N^T (x_,r x x_,s) dr ds
 *
 * in the configuration selected by the follower flag:
 *
 * - dead pressure: reference positions, independent of the displacement,
 * - follower pressure: current positions, additionally contributing the
 *   analytic derivative df/du = df/dx of the quadrature-based force.
 *
 * Because current positions are x = X + u, the derivative with respect to the
 * nodal positions equals the derivative with respect to the translational
 * displacements. A follower pressure without a requested load tangent would be
 * linearized inconsistently and is reported as an error.
 *
 * @param model_data Reference and current nodal geometry and compiled surfaces.
 * @param load Global nodal force field modified in place.
 * @param load_tangent Optional load-tangent triplets receiving df/du.
 * @param system_dof_ids Global equation numbering of the load tangent.
 * @param scale Multiplier of the nominal pressure.
 */
void PLoad::apply(
    model::ModelData&      model_data,
    model::Field&          load,
    TripletList*           load_tangent,
    constraint::Equations&,
    const SystemDofIds&    system_dof_ids,
    TripletList&,
    Precision              scale
) {
    // Validate the definition and the configuration required by the flag
    logging::error(region_ != nullptr,
        "PLOAD: target surface region is not initialized");
    logging::error(std::isfinite(pressure_),
        "PLOAD: pressure must be finite");
    logging::error(model_data.positions != nullptr && model_data.positions_reference != nullptr,
        "PLOAD: positions fields are not initialized");
    logging::error(!follower_ || load_tangent != nullptr,
        "PLOAD: follower pressure requires a load tangent and is only supported in nonlinear analyses");

    const Precision pressure = scale * pressure_;

    if (pressure == Precision(0)) return;

    // Dead pressure acts on the reference surface, follower pressure follows
    // the normal and area of the current surface
    const auto&  positions = follower_ ? *model_data.positions : *model_data.positions_reference;
    TripletList* tangent   = follower_ ? load_tangent : nullptr;

    DynamicVector force;
    DynamicMatrix stiffness;

    for (ID surface_id : *region_) {
        const auto& surface = model_data.surfaces[surface_id];
        if (!surface) continue;

        // Integrate the local nodal force and, for follower pressure, its
        // derivative with respect to the nodal positions
        surface->integrate_pressure(positions, pressure, force, tangent ? &stiffness : nullptr);

        // Scatter the local force into the global nodal load
        const ID* nodes = surface->nodes();

        for (Index i = 0; i < surface->n_nodes; ++i) {
            for (Dim d = 0; d < 3; ++d) {
                load(nodes[i], d) += force(3 * i + d);
            }
        }

        if (!tangent) continue;

        // Scatter the 3x3 nodal blocks of df/dx into the global load tangent
        for (Index i = 0; i < surface->n_nodes; ++i) {
            for (Index k = 0; k < surface->n_nodes; ++k) {
                const Mat3 block = stiffness.block<3, 3>(3 * i, 3 * k);
                add_load_tangent_block(*tangent, system_dof_ids, nodes[i], 0, nodes[k], block);
            }
        }
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
