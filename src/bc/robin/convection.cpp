/**
 * @file convection.cpp
 * @brief Implements linear thermal convection as a Robin boundary condition.
 */

#include "convection.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>

namespace fem::bc {

void Convection::apply(model::ModelData&  model_data,
                       model::Field&      rhs,
                       const SystemDofIds& system_dof_ids,
                       TripletList&        lhs,
                       Precision           time,
                       bool                ignore_amplitude) {
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

    if (h == Precision(0)) return;

    const auto& positions = *model_data.positions_reference;
    const Precision source = h * ambient_temperature_;

    for (ID surface_id : *region_) {
        logging::error(surface_id >= 0
                    && static_cast<Index>(surface_id) < model_data.surfaces.size(),
            "CONVECTION: surface ", surface_id, " is outside the compiled surface domain");

        const auto& surface = model_data.surfaces[static_cast<Index>(surface_id)];
        logging::error(surface != nullptr,
            "CONVECTION: surface ", surface_id, " is not initialized");

        // f_h = integral h*T_inf*N dGamma. This uses the ordinary surface rule.
        surface->integrate_scalar_field(
            positions,
            rhs,
            [source](const Vec3&) -> Precision { return source; }
        );

        // K_h = integral h*N*N^T dGamma. Shape products use their dedicated,
        // higher-order quadrature without changing ordinary surface integration.
        const DynamicMatrix local = surface->integrate_scalar_shape_matrix(
            positions,
            [h](const Vec3&) -> Precision { return h; }
        );

        logging::error(local.rows() == static_cast<Eigen::Index>(surface->n_nodes)
                    && local.cols() == static_cast<Eigen::Index>(surface->n_nodes),
            "CONVECTION: local Robin matrix does not match surface connectivity");
        logging::error(local.allFinite(),
            "CONVECTION: local Robin matrix contains NaN or Inf");

        for (Index i = 0; i < surface->n_nodes; ++i) {
            const ID node_i = surface->nodes()[i];
            const int row = system_dof_ids(static_cast<Eigen::Index>(node_i), 0);
            logging::error(row >= 0,
                "CONVECTION: surface references thermally inactive node ", node_i);

            for (Index j = 0; j < surface->n_nodes; ++j) {
                const ID node_j = surface->nodes()[j];
                const int col = system_dof_ids(static_cast<Eigen::Index>(node_j), 0);
                logging::error(col >= 0,
                    "CONVECTION: surface references thermally inactive node ", node_j);

                const Precision value = local(
                    static_cast<Eigen::Index>(i),
                    static_cast<Eigen::Index>(j)
                );
                if (value != Precision(0)) {
                    lhs.emplace_back(row, col, value);
                }
            }
        }
    }
}

std::string Convection::str() const {
    std::ostringstream os;
    os << "CONVECTION: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", h=" << film_coefficient_
       << ", ambient=" << ambient_temperature_;

    if (amplitude_) os << ", amplitude=" << amplitude_->name;
    return os.str();
}

} // namespace fem::bc
