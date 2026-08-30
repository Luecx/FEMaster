/**
 * @file thermal_surface.h
 * @brief Provides consistent scalar surface integration for thermal Neumann data.
 */

#pragma once

#include "../../core/logging.h"
#include "../../core/types_eig.h"
#include "../../data/field.h"
#include "../../model/geometry/surface/surface3.h"
#include "../../model/geometry/surface/surface4.h"
#include "../../model/geometry/surface/surface6.h"
#include "../../model/geometry/surface/surface8.h"

#include <utility>

namespace fem::bc::detail {

template<class Callback>
void visit_thermal_surface(const model::SurfaceInterface& surface, Callback&& callback) {
    if (const auto* typed = dynamic_cast<const model::Surface3*>(&surface)) {
        callback(*typed);
        return;
    }
    if (const auto* typed = dynamic_cast<const model::Surface4*>(&surface)) {
        callback(*typed);
        return;
    }
    if (const auto* typed = dynamic_cast<const model::Surface6*>(&surface)) {
        callback(*typed);
        return;
    }
    if (const auto* typed = dynamic_cast<const model::Surface8*>(&surface)) {
        callback(*typed);
        return;
    }

    logging::error(false,
        "Thermal surface load: unsupported surface topology with ", surface.n_nodes, " nodes");
}

template<Index N>
void assemble_flux_on_surface(const model::Surface<N>& surface,
                              const model::Field&      positions,
                              const SystemDofIds&      thermal_dof_ids,
                              Precision                flux,
                              DynamicVector&           heat_source) {
    const auto coordinates = surface.node_coords_global(positions);
    const auto& scheme      = surface.integration_scheme();

    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);
        const auto shape = surface.shape_function(point.r, point.s);
        const auto jac   = surface.jacobian(coordinates, point.r, point.s);
        const Precision weighted_area = jac.col(0).cross(jac.col(1)).norm() * point.w;

        for (Index i = 0; i < N; ++i) {
            const ID node = surface.nodeIds[i];
            const int dof = thermal_dof_ids(static_cast<Index>(node), 0);
            logging::error(dof >= 0,
                "Thermal surface load targets thermally inactive node ", node);
            heat_source(dof) += shape(i) * flux * weighted_area;
        }
    }
}

template<Index N>
void assemble_convection_on_surface(const model::Surface<N>& surface,
                                    const model::Field&      positions,
                                    const SystemDofIds&      thermal_dof_ids,
                                    Precision                film,
                                    Precision                ambient_temperature,
                                    TripletList&              matrix_terms,
                                    DynamicVector&           heat_source) {
    const auto coordinates = surface.node_coords_global(positions);
    const auto& scheme      = surface.integration_scheme();

    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);
        const auto shape = surface.shape_function(point.r, point.s);
        const auto jac   = surface.jacobian(coordinates, point.r, point.s);
        const Precision weighted_area = jac.col(0).cross(jac.col(1)).norm() * point.w;

        for (Index i = 0; i < N; ++i) {
            const ID node_i = surface.nodeIds[i];
            const int dof_i = thermal_dof_ids(static_cast<Index>(node_i), 0);
            logging::error(dof_i >= 0,
                "Convection targets thermally inactive node ", node_i);

            heat_source(dof_i) += shape(i) * film * ambient_temperature * weighted_area;

            for (Index j = 0; j < N; ++j) {
                const ID node_j = surface.nodeIds[j];
                const int dof_j = thermal_dof_ids(static_cast<Index>(node_j), 0);
                logging::error(dof_j >= 0,
                    "Convection targets thermally inactive node ", node_j);
                matrix_terms.emplace_back(
                    dof_i,
                    dof_j,
                    shape(i) * shape(j) * film * weighted_area
                );
            }
        }
    }
}

} // namespace fem::bc::detail
