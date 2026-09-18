/**
 * @file convection.cpp
 * @brief Implements linear thermal convection RHS and boundary-operator assembly.
 *
 * Newton cooling separates naturally into a prescribed ambient source
 *
 *     q_h = integral_Gamma h T_inf N^T dGamma
 *
 * and an unknown-dependent symmetric boundary operator
 *
 *     K_h = integral_Gamma h N N^T dGamma.
 *
 * Both terms are integrated in the reference configuration. Shape-product
 * integration uses a dedicated higher-order surface quadrature because
 * `N_i N_j` has higher polynomial order than ordinary surface loading.
 *
 * @see Convection
 * @see model::SurfaceInterface::integrate_scalar_shape_matrix
 *
 * @author Finn Eggers
 * @date 18.09.2026
 */

#include "convection.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <cmath>
#include <cstddef>
#include <sstream>

namespace fem::bc {

/**
 * Evaluates the effective film coefficient used by both convection contributions.
 *
 * The nominal film coefficient must be finite and non-negative. If an amplitude
 * is active, it scales `h`; the resulting effective coefficient is validated
 * again because an arbitrary amplitude could otherwise reverse the dissipative
 * boundary operator or introduce a non-finite value.
 *
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Use the nominal coefficient directly when true.
 * @return Validated effective film coefficient.
 */
Precision Convection::effective_film_coefficient(Precision time, bool ignore_amplitude) const {
    // Validate the physical parameters before applying temporal scaling
    logging::error(std::isfinite(film_coefficient_) && film_coefficient_ >= Precision(0),
        "CONVECTION: film coefficient must be finite and non-negative");
    logging::error(std::isfinite(ambient_temperature_),
        "CONVECTION: ambient temperature must be finite");

    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    const Precision h = film_coefficient_ * scale;

    logging::error(std::isfinite(h) && h >= Precision(0),
        "CONVECTION: effective film coefficient must be finite and non-negative");

    return h;
}

/**
 * Assembles the prescribed ambient source into the scalar thermal RHS.
 *
 * The source term of Newton cooling is
 *
 *     q_h = integral_Gamma h T_inf N^T dGamma.
 *
 * It is independent of the unknown temperature field and can therefore be
 * assembled through the ordinary load-like RHS interface.
 *
 * @param model_data Compiled surface topology and reference geometry.
 * @param rhs Scalar nodal thermal RHS receiving the ambient source.
 * @param time Analysis time used for optional amplitude evaluation.
 * @param ignore_amplitude Apply the nominal film coefficient when true.
 */
void Convection::apply(model::ModelData& model_data,
                       model::Field&     rhs,
                       Precision         time,
                       bool              ignore_amplitude) {
    // Validate the target and scalar thermal assembly context
    logging::error(region_ != nullptr,
        "CONVECTION: target surface region is not set");
    logging::error(model_data.positions_reference != nullptr,
        "CONVECTION: reference positions are not initialized");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components == 1,
        "CONVECTION: target field must be a NODE field with exactly one component");

    // A vanishing film coefficient removes both convection contributions exactly
    const Precision h = effective_film_coefficient(time, ignore_amplitude);
    if (h == Precision(0)) {
        return;
    }

    const Precision source = h * ambient_temperature_;

    // Integrate the ambient heat source consistently over every selected surface
    for (ID surface_id : *region_) {
        logging::error(surface_id >= 0
                    && static_cast<Index>(surface_id) < static_cast<Index>(model_data.surfaces.size()),
            "CONVECTION: surface ", surface_id, " is outside the compiled surface domain");

        const auto& surface = model_data.surfaces[static_cast<std::size_t>(surface_id)];
        logging::error(surface != nullptr,
            "CONVECTION: surface ", surface_id, " is not initialized");

        surface->integrate_scalar_field(
            *model_data.positions_reference,
            rhs,
            [source](const Vec3&) -> Precision { return source; }
        );
    }
}

/**
 * Assembles the temperature-dependent convection boundary operator.
 *
 * For each selected surface the local matrix
 *
 *     K_h^e = integral_Gamma_e h N N^T dGamma
 *
 * is integrated in connectivity ordering. The scalar thermal DOF map converts
 * every surface node to an active system row and column before non-zero entries
 * are appended to the global sparse triplet list.
 *
 * @param model_data Compiled surface topology and reference geometry.
 * @param system_dof_ids Scalar node-to-active-temperature equation mapping.
 * @param matrix Sparse triplet list receiving convection operator entries.
 * @param time Analysis time used for optional amplitude evaluation.
 * @param ignore_amplitude Apply the nominal film coefficient when true.
 */
void Convection::apply_matrix(model::ModelData&   model_data,
                              const SystemDofIds& system_dof_ids,
                              TripletList&        matrix,
                              Precision           time,
                              bool                ignore_amplitude) {
    // Validate the surface target and scalar thermal system numbering
    logging::error(region_ != nullptr,
        "CONVECTION: target surface region is not set");
    logging::error(model_data.positions_reference != nullptr,
        "CONVECTION: reference positions are not initialized");
    logging::error(system_dof_ids.rows() == model_data.positions_reference->rows,
        "CONVECTION: thermal DOF map does not match the nodal domain");
    logging::error(system_dof_ids.cols() == 1,
        "CONVECTION: thermal DOF map must contain exactly one component");

    const Precision h = effective_film_coefficient(time, ignore_amplitude);
    if (h == Precision(0)) {
        return;
    }

    const auto& positions = *model_data.positions_reference;

    // Assemble one consistent Robin matrix for every selected reference surface
    for (ID surface_id : *region_) {
        logging::error(surface_id >= 0
                    && static_cast<Index>(surface_id) < static_cast<Index>(model_data.surfaces.size()),
            "CONVECTION: surface ", surface_id, " is outside the compiled surface domain");

        const auto& surface = model_data.surfaces[static_cast<std::size_t>(surface_id)];
        logging::error(surface != nullptr,
            "CONVECTION: surface ", surface_id, " is not initialized");

        // Higher-order shape-product quadrature evaluates integral h N N^T dGamma
        const DynamicMatrix local = surface->integrate_scalar_shape_matrix(
            positions,
            [h](const Vec3&) -> Precision { return h; }
        );

        logging::error(local.rows() == surface->n_nodes && local.cols() == surface->n_nodes,
            "CONVECTION: local boundary matrix does not match surface connectivity");
        logging::error(local.allFinite(),
            "CONVECTION: local boundary matrix contains NaN or Inf");

        // Map the local surface operator to active scalar thermal system indices
        for (Index i = 0; i < surface->n_nodes; ++i) {
            const ID  node_i = surface->nodes()[i];
            const int row    = system_dof_ids(static_cast<Eigen::Index>(node_i), 0);

            logging::error(row >= 0,
                "CONVECTION: surface references thermally inactive node ", node_i);

            for (Index j = 0; j < surface->n_nodes; ++j) {
                const ID  node_j = surface->nodes()[j];
                const int col    = system_dof_ids(static_cast<Eigen::Index>(node_j), 0);

                logging::error(col >= 0,
                    "CONVECTION: surface references thermally inactive node ", node_j);

                const Precision value = local(
                    static_cast<Eigen::Index>(i),
                    static_cast<Eigen::Index>(j)
                );

                // Adjacent surfaces may contribute to the same entry. Duplicate
                // triplets are intentionally left for the global sparse assembly
                // to sum.
                if (value != Precision(0)) {
                    matrix.emplace_back(row, col, value);
                }
            }
        }
    }
}

/**
 * Builds the diagnostic representation of the convection condition.
 *
 * @return Human-readable target region, nominal film coefficient, ambient
 *         temperature and optional amplitude.
 */
std::string Convection::str() const {
    std::ostringstream os;

    os << "CONVECTION: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", h=" << film_coefficient_
       << ", ambient=" << ambient_temperature_;

    if (amplitude_) {
        os << ", amplitude=" << amplitude_->name;
    }

    return os.str();
}

} // namespace fem::bc
