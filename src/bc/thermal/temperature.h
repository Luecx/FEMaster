/**
 * @file temperature.h
 * @brief Defines prescribed-temperature thermal conditions.
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
 * @see Condition
 * @see constraint::Equation
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include "../condition.h"
#include "../../constraints/types/equation.h"
#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../data/field.h"
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
 * The class defines only the prescribed thermal primary variable. Heat flux and
 * convection are independent thermal condition classes with their own
 * solver-facing operations.
 *
 * Condition-history identity is the semantic target region. The prescribed
 * scalar temperature is the value of that condition and may be replaced without
 * changing its identity.
 */
struct Temperature : Condition {
    // Types
    using Ptr = std::shared_ptr<Temperature>;

    // Exactly one target region defines the nodes receiving the temperature
    model::NodeRegion::Ptr    node_region_    = nullptr;
    model::SurfaceRegion::Ptr surface_region_ = nullptr;
    model::ElementRegion::Ptr element_region_ = nullptr;

    // Absolute prescribed temperature defining the scalar equation T_i = T_bar
    Precision temperature_ = Precision(0);

    // Construction from node, surface or element targets
    Temperature(model::NodeRegion::Ptr node_region, Precision temperature)
        : node_region_(std::move(node_region)), temperature_(temperature) {}
    Temperature(model::SurfaceRegion::Ptr surface_region, Precision temperature)
        : surface_region_(std::move(surface_region)), temperature_(temperature) {}
    Temperature(model::ElementRegion::Ptr element_region, Precision temperature)
        : element_region_(std::move(element_region)), temperature_(temperature) {}
    ~Temperature() override = default;

    // Resolve the configured node, surface or element region to its global nodes
    // and append exactly one prescribed-temperature equation for each unique node.
    //
    // For every selected node i, the scalar thermal DOF receives the equation
    //
    //     1 * T_i = temperature_
    //
    // Surface and element regions are expanded through their connectivity;
    // shared nodes are deduplicated before equations are appended.
    //
    // This condition prescribes an absolute temperature, not a thermal load.
    // It does not interpolate temperature histories or evaluate amplitudes:
    // time, ignore_amplitude and step_progress are unused.
    //
    // Only equations is modified; rhs, system_dof_ids and matrix are unchanged.
    void apply(
        model::ModelData&      model_data,
        model::Field&          rhs,
        constraint::Equations& equations,
        const SystemDofIds&    system_dof_ids,
        TripletList&           matrix,
        Precision              time,
        bool                   ignore_amplitude = false,
        Precision              step_progress   = Precision(1)
    ) override;

    // Diagnostics
    std::string str() const override;
};

} // namespace fem::bc
