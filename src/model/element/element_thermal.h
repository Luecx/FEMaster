/**
 * @file element_thermal.h
 * @brief Declares the thermal finite-element capability interface.
 *
 * Thermal behavior is modeled as an element capability alongside structural
 * mechanics. Concrete elements may implement both interfaces; thermal model
 * assembly selects `ThermalElement` implementations when conductivity, capacity
 * or conductive heat-flux operations are requested.
 *
 * Steady-state analysis consumes conductivity and heat flux. Capacity remains
 * part of the interface for a future transient thermal formulation.
 *
 * @see ThermalElement
 * @see SolidElement
 * @see loadcase::SteadyStateThermal
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "element.h"

namespace fem::model {

/**
 * @brief Provides scalar heat-transfer operations for a finite element.
 *
 * A thermal element associates one scalar temperature unknown with every
 * connectivity node. `conductivity()` and `capacity()` return element matrices
 * in that nodal order using caller-owned storage. `compute_heat_flux()`
 * differentiates a converged nodal temperature field at the element's natural
 * node coordinates and writes the resulting global Cartesian flux into the
 * compiled element-nodal rows.
 */
struct ThermalElement {
    // Polymorphic destruction through the thermal capability interface
    virtual ~ThermalElement() = default;

    // Scalar nodal conductivity and capacity matrices. Returned maps reference
    // the supplied caller-owned storage.
    virtual MapMatrix conductivity(Precision* buffer) = 0;
    virtual MapMatrix capacity    (Precision* buffer) = 0;

    // Recover global conductive heat flux in the compiled element-nodal range
    virtual void compute_heat_flux(Field& heat_flux, const Field& temperature) = 0;
};

} // namespace fem::model
