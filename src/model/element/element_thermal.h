
/**
 * @file element_thermal.h
 * @brief Declares the thermal finite-element capability interface.
 *
 * Thermal behavior is modeled as an orthogonal element capability alongside
 * structural mechanics. Concrete solid elements may implement both interfaces;
 * model assembly cross-casts the common `ElementInterface` pointer to
 * `ThermalElement` when conductivity, capacity or heat-flux operations are
 * requested.
 *
 * Steady-state analysis consumes only conductivity and heat flux. Capacity is
 * retained for a future transient thermal load case.
 *
 * @see ThermalElement
 * @see SolidElement
 * @see loadcase::SteadyStateThermal
 *
 * @author Finn Eggers
 * @date 29.08.2026
 */

#pragma once

#include "element.h"

namespace fem::model {

/**
 * @brief Provides scalar heat-transfer operations for a finite element.
 *
 * A thermal element associates one scalar temperature unknown with every
 * connectivity node. `conductivity()` and `capacity()` return element matrices
 * in that nodal ordering using caller-owned storage. `compute_heat_flux()`
 * differentiates a converged global nodal temperature field and writes global
 * Cartesian flux components to the element's compiled integration-point range.
 *
 * Implementations must not own the supplied matrix buffer or result fields.
 * Their topology, section, material and compiled offsets remain owned by the
 * accompanying `ElementInterface` subobject of the concrete element.
 */
struct ThermalElement {
    // Polymorphic destruction through the thermal capability interface
    virtual ~ThermalElement() = default;

    // Scalar nodal conductivity and capacity matrices in element connectivity
    // order. The returned maps reference the caller-owned buffer.
    virtual MapMatrix conductivity(Precision* buffer) = 0;
    virtual MapMatrix capacity    (Precision* buffer) = 0;

    // Recover global Cartesian heat flux into the element's compiled IP rows
    virtual void compute_heat_flux(Field& heat_flux, const Field& temperature) = 0;
};

} // namespace fem::model
