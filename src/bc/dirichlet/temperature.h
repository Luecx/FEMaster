/**
 * @file temperature.h
 * @brief Defines prescribed-temperature Dirichlet conditions.
 *
 * A temperature condition assigns one scalar thermal primary variable to nodes
 * selected directly or through element and surface connectivity. Application
 * converts the semantic target into nodal equations of the form
 *
 *     T_i = T_bar.
 *
 * Shared nodes from element or surface targets are deduplicated before the
 * equations are emitted so one physical node receives one scalar prescription
 * from a single `Temperature` object.
 *
 * @see Dirichlet
 * @see ThermalCondition
 * @see constraint::Equation
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include "dirichlet.h"
#include "../thermal.h"
#include "../../core/types_num.h"
#include "../../data/region.h"

#include <memory>
#include <string>
#include <utility>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Prescribes one absolute temperature on a model region.
 *
 * Exactly one of `node_region_`, `surface_region_` or `element_region_` defines
 * the semantic target. Surface and element targets are expanded to their global
 * connected nodes during application.
 *
 * FEMaster represents the scalar thermal unknown as local thermal DOF zero. For
 * every unique resolved node `i`, the condition therefore appends
 *
 *     1 * T_i = temperature_.
 *
 * The class defines only the prescribed primary variable. Thermal fluxes and
 * convection belong to the Neumann and mixed categories respectively.
 */
struct Temperature : Dirichlet, ThermalCondition {
    // Types
    using Ptr = std::shared_ptr<Temperature>;

    // Exactly one target region defines the nodes receiving the temperature
    model::NodeRegion::Ptr    node_region_    = nullptr;
    model::SurfaceRegion::Ptr surface_region_ = nullptr;
    model::ElementRegion::Ptr element_region_ = nullptr;

    // Absolute prescribed temperature used as the scalar equation RHS
    Precision temperature_ = Precision(0);

    // Construction from node, surface or element targets
    Temperature(model::NodeRegion::Ptr node_region, Precision temperature)
        : node_region_(std::move(node_region)), temperature_(temperature) {}
    Temperature(model::SurfaceRegion::Ptr surface_region, Precision temperature)
        : surface_region_(std::move(surface_region)), temperature_(temperature) {}
    Temperature(model::ElementRegion::Ptr element_region, Precision temperature)
        : element_region_(std::move(element_region)), temperature_(temperature) {}
    ~Temperature() override = default;

    // Resolve the target nodes, remove duplicates and append T_i = T_bar rows
    void apply(model::ModelData& model_data, constraint::Equations& equations) override;

    // Diagnostics
    std::string str() const override;
};

} // namespace fem::bc
