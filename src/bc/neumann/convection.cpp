/**
 * @file convection.cpp
 * @brief Implements thermal convection boundary conditions.
 *
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

void Convection::apply(model::ModelData& model_data,
                       model::Field&     bc,
                       Precision         time,
                       bool              ignore_amplitude) {
    logging::error(region_ != nullptr,
        "CONVECTION: target surface region is not set");
    logging::error(model_data.positions_reference != nullptr,
        "CONVECTION: reference positions are not initialized");
    logging::error(bc.domain == model::FieldDomain::NODE && bc.components >= 3,
        "CONVECTION: target field must be a NODE field with at least three components");
    logging::error(std::isfinite(film_coefficient_) && film_coefficient_ >= Precision(0),
        "CONVECTION: film coefficient must be finite and non-negative");
    logging::error(std::isfinite(ambient_temperature_),
        "CONVECTION: ambient temperature must be finite");

    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    const Vec3 value{
        film_coefficient_ * ambient_temperature_ * scale,
        Precision(0),
        Precision(0)
    };

    for (ID surface_id : *region_) {
        logging::error(surface_id >= 0
                    && static_cast<Index>(surface_id) < model_data.surfaces.size(),
            "CONVECTION: surface ", surface_id, " is outside the compiled surface domain");

        const auto& surface = model_data.surfaces[static_cast<Index>(surface_id)];
        logging::error(surface != nullptr,
            "CONVECTION: surface ", surface_id, " is not initialized");

        surface->integrate_vector_field(
            *model_data.positions_reference,
            bc,
            [value](const Vec3&) -> Vec3 { return value; }
        );
    }
}

void Convection::apply_conductivity(model::ModelData& model_data,
                                    const SystemDofIds& thermal_dof_ids,
                                    TripletList&         triplets,
                                    Precision            time,
                                    bool                 ignore_amplitude) {
    logging::error(region_ != nullptr,
        "CONVECTION: target surface region is not set");
    logging::error(model_data.positions_reference != nullptr,
        "CONVECTION: reference positions are not initialized");
    logging::error(thermal_dof_ids.rows() == model_data.positions_reference->rows,
        "CONVECTION: thermal DOF map does not match the nodal domain");
    logging::error(std::isfinite(film_coefficient_) && film_coefficient_ >= Precision(0),
        "CONVECTION: film coefficient must be finite and non-negative");

    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    const Precision h = film_coefficient_ * scale;

    if (h == Precision(0))
        return;

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
                    const int dof_i = thermal_dof_ids(static_cast<Index>(node_i), 0);
                    if (dof_i < 0)
                        continue;

                    for (Index j = 0; j < surface->n_nodes; ++j) {
                        const ID node_j = surface->nodes()[j];
                        const int dof_j = thermal_dof_ids(static_cast<Index>(node_j), 0);
                        if (dof_j < 0)
                            continue;

                        triplets.emplace_back(
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
