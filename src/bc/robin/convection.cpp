/**
 * @file convection.cpp
 * @brief Implements linear thermal convection as a Robin boundary condition.
 *
 * Convection contributes a prescribed ambient source and a temperature-dependent
 * boundary operator. The source is integrated with the ordinary surface rule,
 * while the `N N^T` operator term uses the dedicated scalar shape-product
 * integration routine so its higher polynomial order does not change the
 * established quadrature used by ordinary surface loads.
 *
 * Both terms are evaluated in the reference configuration. The local Robin
 * matrix is mapped directly to active scalar temperature DOFs and appended as
 * sparse triplets for later global assembly.
 *
 * @see convection.h
 * @see model::SurfaceInterface::integrate_scalar_shape_matrix
 *
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
 * Assembles the consistent ambient source and Robin boundary operator.
 *
 * The effective film coefficient is obtained from the nominal coefficient and
 * optional amplitude. For each selected surface, the ambient term
 * `h T_inf` is integrated consistently into the scalar nodal RHS, and the local
 * matrix `integral h N N^T dGamma` is evaluated with shape-product quadrature.
 * The local matrix entries are then mapped through the scalar thermal DOF matrix
 * and appended to the global triplet list.
 *
 * A zero effective film coefficient makes both contributions vanish and returns
 * before any surface traversal. Every surface node must be thermally active;
 * otherwise a valid local Robin matrix could not be represented in the supplied
 * algebraic thermal system.
 *
 * @param model_data Compiled surfaces and reference nodal geometry.
 * @param rhs Scalar nodal thermal right-hand side receiving `h T_inf` terms.
 * @param system_dof_ids Scalar node-to-system temperature DOF mapping.
 * @param lhs Sparse global triplet list receiving the Robin operator terms.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether the nominal film coefficient is used unscaled.
 */
void Convection::apply(model::ModelData&   model_data,
                       model::Field&       rhs,
                       const SystemDofIds& system_dof_ids,
                       TripletList&         lhs,
                       Precision           time,
                       bool                ignore_amplitude) {
    // Validate the target region, reference geometry and scalar thermal assembly
    // layouts before evaluating any physical convection parameters.
    logging::error(region_ != nullptr,
        "CONVECTION: target surface region is not set");
    logging::error(model_data.positions_reference != nullptr,
        "CONVECTION: reference positions are not initialized");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components == 1,
        "CONVECTION: target field must be a NODE field with exactly one component");
    logging::error(system_dof_ids.rows() == static_cast<Eigen::Index>(model_data.positions_reference->rows),
        "CONVECTION: thermal DOF map does not match the nodal domain");
    logging::error(system_dof_ids.cols() == 1,
        "CONVECTION: thermal DOF map must contain exactly one component");

    // Validate the nominal physical parameters before optional temporal scaling.
    // Negative film coefficients would reverse the dissipative Robin operator and
    // are not supported by the linear convection model.
    logging::error(std::isfinite(film_coefficient_) && film_coefficient_ >= Precision(0),
        "CONVECTION: film coefficient must be finite and non-negative");
    logging::error(std::isfinite(ambient_temperature_),
        "CONVECTION: ambient temperature must be finite");

    // Apply the optional amplitude to the film coefficient. The ambient
    // temperature remains the prescribed environmental state and is not scaled.
    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    const Precision h = film_coefficient_ * scale;

    // Revalidate the effective coefficient because a user-defined amplitude can
    // itself generate non-finite or negative values even when the nominal value
    // is valid.
    logging::error(std::isfinite(h) && h >= Precision(0),
        "CONVECTION: effective film coefficient must be finite and non-negative");

    // With h = 0 both the ambient source and temperature-dependent boundary
    // operator vanish exactly.
    if (h == Precision(0)) return;

    // Convection is integrated on the reference boundary. The scalar source
    // `h T_inf` is constant over every surface for the current implementation.
    const auto& positions = *model_data.positions_reference;
    const Precision source = h * ambient_temperature_;

    // Assemble both Robin contributions independently for every selected compiled
    // surface.
    for (ID surface_id : *region_) {
        // Validate and resolve the dense compiled surface identifier before
        // invoking geometric integration routines.
        logging::error(surface_id >= 0
                    && static_cast<Index>(surface_id) < model_data.surfaces.size(),
            "CONVECTION: surface ", surface_id, " is outside the compiled surface domain");

        const auto& surface = model_data.surfaces[static_cast<Index>(surface_id)];
        logging::error(surface != nullptr,
            "CONVECTION: surface ", surface_id, " is not initialized");

        // Assemble f_h = integral_Gamma h T_inf N dGamma into the scalar nodal
        // thermal RHS. This linear-in-N integrand uses the ordinary surface rule.
        surface->integrate_scalar_field(
            positions,
            rhs,
            [source](const Vec3&) -> Precision { return source; }
        );

        // Assemble the local Robin operator K_h = integral_Gamma h N N^T dGamma.
        // Products of interpolation functions have a higher polynomial order, so
        // this operation deliberately uses its dedicated shape-matrix quadrature.
        const DynamicMatrix local = surface->integrate_scalar_shape_matrix(
            positions,
            [h](const Vec3&) -> Precision { return h; }
        );

        // The local operator must match the surface connectivity exactly and must
        // remain finite before any entries are added to the global sparse system.
        logging::error(local.rows() == static_cast<Eigen::Index>(surface->n_nodes)
                    && local.cols() == static_cast<Eigen::Index>(surface->n_nodes),
            "CONVECTION: local Robin matrix does not match surface connectivity");
        logging::error(local.allFinite(),
            "CONVECTION: local Robin matrix contains NaN or Inf");

        // Map every local row node to its active scalar thermal system DOF.
        for (Index i = 0; i < surface->n_nodes; ++i) {
            const ID node_i = surface->nodes()[i];
            const int row   = system_dof_ids(static_cast<Eigen::Index>(node_i), 0);
            logging::error(row >= 0,
                "CONVECTION: surface references thermally inactive node ", node_i);

            // Map every local column node independently and append non-zero
            // entries. Duplicate triplets from adjacent surfaces are intentionally
            // left for the global sparse constructor to sum.
            for (Index j = 0; j < surface->n_nodes; ++j) {
                const ID node_j = surface->nodes()[j];
                const int col   = system_dof_ids(static_cast<Eigen::Index>(node_j), 0);
                logging::error(col >= 0,
                    "CONVECTION: surface references thermally inactive node ", node_j);

                const Precision value = local(
                    static_cast<Eigen::Index>(i),
                    static_cast<Eigen::Index>(j)
                );

                // Skip exact zeros to keep the sparse triplet list compact without
                // applying an arbitrary numerical threshold to a physical matrix.
                if (value != Precision(0)) {
                    lhs.emplace_back(row, col, value);
                }
            }
        }
    }
}

/**
 * Builds the diagnostic representation of the convection condition.
 *
 * The result reports the target surface region, nominal film coefficient and
 * ambient temperature before optional amplitude scaling.
 *
 * @return Human-readable convection description.
 */
std::string Convection::str() const {
    std::ostringstream os;

    // Report the region identity and the two physical convection parameters.
    os << "CONVECTION: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", h=" << film_coefficient_
       << ", ambient=" << ambient_temperature_;

    // Include temporal scaling only when an amplitude is assigned.
    if (amplitude_) os << ", amplitude=" << amplitude_->name;

    return os.str();
}

} // namespace fem::bc
