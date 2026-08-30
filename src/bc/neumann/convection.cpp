/**
 * @file convection.cpp
 * @brief Implements linear thermal convection boundary conditions.
 *
 * Convection is assembled as one load-side operation even though it contributes
 * to both sides of the linear thermal system. The ambient-temperature part is
 * integrated into component zero of the nodal RHS using the existing surface
 * vector integrator. The temperature-dependent film term is integrated into the
 * sparse LHS through the same `apply()` call.
 *
 * @see convection.h
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "convection.h"

#include "../../core/logging.h"
#include "../../math/quadrature.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>

namespace fem::bc {

/**
 * Assembles the complete linear convection boundary contribution.
 *
 * Starting from \f$q = h(T_\infty - T)\f$, finite-element discretization gives
 * the ambient source term
 * \f$\mathbf{f}_h = \int_\Gamma hT_\infty\mathbf{N}^T\,\mathrm{d}\Gamma\f$
 * and the boundary conductivity term
 * \f$\mathbf{K}_h = \int_\Gamma h\mathbf{N}^T\mathbf{N}\,\mathrm{d}\Gamma\f$.
 * Both use the same amplitude-scaled film coefficient and reference geometry.
 *
 * The RHS is represented by FEMaster's existing three-component nodal load field
 * and only component zero is used for thermal power. The LHS is assembled using
 * scalar thermal DOF identifiers from column zero of `system_dof_ids`.
 *
 * @param model_data Compiled surface topology and reference nodal positions.
 * @param rhs Three-component nodal thermal RHS. Only component zero is modified.
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
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components >= 3,
        "CONVECTION: target field must be a NODE field with at least three components");
    logging::error(system_dof_ids != nullptr,
        "CONVECTION: system DOF map is required for LHS assembly");
    logging::error(lhs != nullptr,
        "CONVECTION: LHS triplet list is required for boundary-matrix assembly");
    logging::error(system_dof_ids->rows() == model_data.positions_reference->rows,
        "CONVECTION: thermal DOF map does not match the nodal domain");
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

    // The existing surface-vector integration is reused for the scalar thermal
    // source by storing h*T_inf exclusively in component zero.
    const Vec3 rhs_value{
        h * ambient_temperature_,
        Precision(0),
        Precision(0)
    };

    // Triangular subdivision handles both triangular and quadrilateral surface
    // polygons through one common surface-integration path.
    const math::quadrature::Quadrature scheme(
        math::quadrature::DOMAIN_ISO_TRI,
        math::quadrature::ORDER_CUBIC
    );

    for (ID surface_id : *region_) {
        logging::error(surface_id >= 0
                    && static_cast<Index>(surface_id) < model_data.surfaces.size(),
            "CONVECTION: surface ", surface_id, " is outside the compiled surface domain");

        const auto& surface = model_data.surfaces[static_cast<Index>(surface_id)];
        logging::error(surface != nullptr,
            "CONVECTION: surface ", surface_id, " is not initialized");

        // Assemble the ambient-temperature source contribution f_h.
        surface->integrate_vector_field(
            *model_data.positions_reference,
            rhs,
            [rhs_value](const Vec3&) -> Vec3 { return rhs_value; }
        );

        // Assemble the temperature-dependent boundary matrix K_h.
        const auto polygon = surface->local_domain_polygon();
        surface->integrate_triangular(
            *model_data.positions_reference,
            polygon,
            scheme,
            [&](const Vec2& local, const Vec3&, Precision weight) {
                const DynamicVector shape = surface->shape_function(local);
                logging::error(shape.size() == static_cast<Eigen::Index>(surface->n_nodes),
                    "CONVECTION: surface shape-function size does not match connectivity");

                for (Index i = 0; i < surface->n_nodes; ++i) {
                    const ID node_i = surface->nodes()[i];
                    const int dof_i = (*system_dof_ids)(static_cast<Index>(node_i), 0);
                    if (dof_i < 0)
                        continue;

                    for (Index j = 0; j < surface->n_nodes; ++j) {
                        const ID node_j = surface->nodes()[j];
                        const int dof_j = (*system_dof_ids)(static_cast<Index>(node_j), 0);
                        if (dof_j < 0)
                            continue;

                        lhs->emplace_back(
                            dof_i,
                            dof_j,
                            h * shape(static_cast<Eigen::Index>(i))
                              * shape(static_cast<Eigen::Index>(j))
                              * weight
                        );
                    }
                }
            }
        );
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
