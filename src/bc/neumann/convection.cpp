/**
 * @file convection.cpp
 * @brief Implements linear thermal convection boundary conditions.
 *
 * Convection is assembled as one load-side operation even though it contributes
 * to both sides of the linear thermal system. The ambient-temperature part is
 * integrated into a scalar nodal RHS. The temperature-dependent film term is
 * integrated over the complete selected surface using the surface rule selected
 * for products of shape functions.
 *
 * @see convection.h
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "convection.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>

namespace fem::bc {

/**
 * Assembles the complete linear convection boundary contribution.
 *
 * Starting from
 *
 *     q = h (T_inf - T),
 *
 * finite-element discretization gives the ambient source term
 *
 *     f_h = integral_Gamma h T_inf N^T dGamma
 *
 * and the boundary conductivity term
 *
 *     K_h = integral_Gamma h N^T N dGamma.
 *
 * Both use the same amplitude-scaled film coefficient and reference geometry.
 *
 * The RHS is a one-component nodal field representing thermal power. Each entry
 * of the symmetric boundary matrix is evaluated with
 * `SurfaceInterface::integrate_scalar_field()`, which integrates over the
 * complete surface using the formulation's full-surface quadrature rule.
 * `integrate_triangular()` is intentionally not used because that routine is for
 * polygon-restricted partial-surface integration.
 *
 * @param model_data Compiled surface topology and reference nodal positions.
 * @param rhs One-component nodal thermal right-hand side.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether amplitude scaling is disabled.
 * @param system_dof_ids Thermal node-to-system DOF map. Must be non-null.
 * @param lhs Sparse triplet list receiving convection boundary-matrix entries.
 *            Must be non-null.
 */
void Convection::apply(model::ModelData&       model_data,
                       model::Field&           rhs,
                       Precision               time,
                       bool                    ignore_amplitude,
                       const SystemDofIds*      system_dof_ids,
                       TripletList*             lhs) {
    logging::error(region_ != nullptr,
        "CONVECTION: target surface region is not set");
    logging::error(model_data.positions_reference != nullptr,
        "CONVECTION: reference positions are not initialized");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components == 1,
        "CONVECTION: target field must be a NODE field with exactly one component");
    logging::error(system_dof_ids != nullptr,
        "CONVECTION: system DOF map is required for LHS assembly");
    logging::error(lhs != nullptr,
        "CONVECTION: LHS triplet list is required for boundary-matrix assembly");
    logging::error(system_dof_ids->rows() == static_cast<Eigen::Index>(model_data.positions_reference->rows),
        "CONVECTION: thermal DOF map does not match the nodal domain");
    logging::error(system_dof_ids->cols() >= 1,
        "CONVECTION: thermal DOF map must contain scalar component zero");
    logging::error(std::isfinite(film_coefficient_) && film_coefficient_ >= Precision(0),
        "CONVECTION: film coefficient must be finite and non-negative");
    logging::error(std::isfinite(ambient_temperature_),
        "CONVECTION: ambient temperature must be finite");

    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    const Precision h = film_coefficient_ * scale;

    if (h == Precision(0))
        return;

    const auto& positions = *model_data.positions_reference;

    const Precision rhs_value = h * ambient_temperature_;

    for (ID surface_id : *region_) {
        logging::error(surface_id >= 0
                    && static_cast<Index>(surface_id) < model_data.surfaces.size(),
            "CONVECTION: surface ", surface_id, " is outside the compiled surface domain");

        const auto& surface = model_data.surfaces[static_cast<Index>(surface_id)];
        logging::error(surface != nullptr,
            "CONVECTION: surface ", surface_id, " is not initialized");

        // Assemble f_h over the complete surface with its ordinary quadrature.
        surface->integrate_scalar_field(
            positions,
            rhs,
            [rhs_value](const Vec3&) -> Precision { return rhs_value; }
        );

        // Assemble only the upper triangle of K_h and mirror off-diagonal terms.
        // The scalar surface integrator owns the full-surface quadrature and
        // physical Jacobian; the natural point is recovered only to evaluate N.
        for (Index i = 0; i < surface->n_nodes; ++i) {
            const ID  node_i = surface->nodes()[i];
            const int dof_i  = (*system_dof_ids)(static_cast<Index>(node_i), 0);
            if (dof_i < 0)
                continue;

            for (Index j = i; j < surface->n_nodes; ++j) {
                const ID  node_j = surface->nodes()[j];
                const int dof_j  = (*system_dof_ids)(static_cast<Index>(node_j), 0);
                if (dof_j < 0)
                    continue;

                const Precision value = surface->integrate_scalar_field(
                    positions,
                    [&](const Vec3& position) -> Precision {
                        const Vec2 local = surface->global_to_local(position, positions);
                        const DynamicVector shape = surface->shape_function(local);

                        logging::error(shape.size() == static_cast<Eigen::Index>(surface->n_nodes),
                            "CONVECTION: surface shape-function size does not match connectivity");

                        return h
                             * shape(static_cast<Eigen::Index>(i))
                             * shape(static_cast<Eigen::Index>(j));
                    }
                );

                lhs->emplace_back(dof_i, dof_j, value);
                if (i != j)
                    lhs->emplace_back(dof_j, dof_i, value);
            }
        }
    }
}

/**
 * Builds the diagnostic representation of the convection boundary condition.
 *
 * @return Human-readable target, film coefficient and ambient temperature.
 */
std::string Convection::str() const {
    std::ostringstream os;
    os << "CONVECTION: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", h=" << film_coefficient_
       << ", ambient=" << ambient_temperature_;

    if (amplitude_)
        os << ", amplitude=" << amplitude_->name;

    return os.str();
}

} // namespace fem::bc
