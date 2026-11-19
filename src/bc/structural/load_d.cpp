/**
 * @file load_d.cpp
 * @brief Implements consistent integration of distributed surface tractions.
 *
 * For every selected surface, the continuous traction field is converted to
 * equivalent nodal forces through the standard weak-form contribution
 *
 *     f_e = integral_Gamma_e N^T t dGamma.
 *
 * The surface implementation owns interpolation, quadrature, the physical
 * area Jacobian and its linearization. `DLoad` provides the traction integrand,
 * including optional scaling and local-to-global coordinate transformation.
 * Curvilinear coordinate systems are evaluated independently at every
 * integration point.
 *
 * Dead tractions act per unit reference area. Follower tractions act per unit
 * current area and contribute the analytic load tangent df/du.
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
 * The calling analysis supplies the multiplier; omitted NaN components
 * contribute zero. Each surface evaluates the consistent nodal force
 *
 *     f_e = integral_Gamma N^T t(x) dGamma
 *
 * in the configuration selected by the follower flag:
 *
 * - dead traction: force per unit reference area, evaluated on the reference
 *   positions and independent of the displacement,
 * - follower traction: force per unit current area, evaluated on the current
 *   positions. The analytic load tangent contains the change of the area
 *   measure and, for a position-dependent orientation, of the direction
 *
 *       dt/dx = d(A(x) t_local)/dx.
 *
 * A follower traction without a requested load tangent would be linearized
 * inconsistently and is reported as an error.
 *
 * @param model_data Reference and current nodal geometry and compiled surfaces.
 * @param load Global nodal force field modified in place.
 * @param load_tangent Optional load-tangent triplets receiving df/du.
 * @param system_dof_ids Global equation numbering of the load tangent.
 * @param scale Multiplier of the nominal traction.
 */
void DLoad::apply(
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
        "DLOAD: target surface region is not initialized");
    logging::error(model_data.positions != nullptr && model_data.positions_reference != nullptr,
        "DLOAD: positions fields are not initialized");
    logging::error(!follower_ || load_tangent != nullptr,
        "DLOAD: follower traction requires a load tangent and is only supported in nonlinear analyses");

    // Omitted components contribute no traction
    Vec3 traction = Vec3::Zero();
    for (Dim i = 0; i < 3; ++i) {
        if (std::isfinite(values_[i])) traction[i] = scale * values_[i];
    }

    if (traction.isZero()) return;

    // Dead traction acts on the reference surface, follower traction on the
    // current surface
    const auto&  positions = follower_ ? *model_data.positions : *model_data.positions_reference;
    TripletList* tangent   = follower_ ? load_tangent : nullptr;

    // Global traction and its spatial gradient. A local traction is rotated
    // into the basis evaluated at every quadrature point because curvilinear
    // axes vary across the surface.
    const VecField traction_field = [&](const Vec3& position) -> Vec3 {
        if (!orientation_) return traction;

        const Vec3 local_point = orientation_->to_local(position);
        const auto axes        = orientation_->get_axes(local_point);

        return axes * traction;
    };
    const TenField traction_gradient = [&](const Vec3& position) -> Mat3 {
        return orientation_->get_axes_gradient(position, traction);
    };

    DynamicVector force;
    DynamicMatrix stiffness;

    for (ID surface_id : *region_) {
        const auto& surface = model_data.surfaces[surface_id];
        if (!surface) continue;

        // Integrate the local nodal force and, for follower traction, its
        // derivative with respect to the nodal positions
        surface->integrate_traction(
            positions,
            traction_field,
            orientation_ ? &traction_gradient : nullptr,
            force,
            tangent ? &stiffness : nullptr
        );

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
