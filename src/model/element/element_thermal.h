/**
 * @file element_thermal.h
 * @brief Declares the thermal finite-element capability interface.
 *
 * Thermal elements provide the local operators required by heat-transfer
 * analyses and recover conductive heat flux from a scalar nodal temperature
 * field. The interface is intentionally independent of structural element
 * behavior so one concrete element may participate in both mechanical and
 * thermal analyses.
 *
 * Conductivity and capacity matrices are element-local scalar operators with one
 * temperature DOF per node. Heat-flux recovery writes one vector value per
 * element node into globally enumerated `ELEMENT_NODAL` storage. The Model
 * layer subsequently projects those discontinuous element-local nodal values to
 * unique global nodes using the same averaging path as structural nodal results.
 *
 * @see ElementInterface
 * @see Model::compute_heat_flux
 *
 * @author Finn Eggers
 * @date 18.09.2026
 */

#pragma once

#include "element.h"

namespace fem::model {

/**
 * @brief Capability implemented by elements participating in thermal analyses.
 *
 * The local conductivity operator represents
 *
 *     K_T^e = integral_Omega grad(N)^T k grad(N) dOmega,
 *
 * while the capacity operator represents
 *
 *     C_T^e = integral_Omega rho c_p N^T N dOmega.
 *
 * Heat-flux recovery evaluates Fourier's law
 *
 *     q = -k grad(T)
 *
 * internally at formulation-safe recovery points and stores the recovered
 * values in the element's disjoint `ELEMENT_NODAL` row range. No global
 * integration-point heat-flux field is required for nodal output.
 */
struct ThermalElement {
    virtual ~ThermalElement() = default;

    // Element-local scalar thermal operators
    virtual MapMatrix conductivity(Precision* buffer) = 0;
    virtual MapMatrix capacity    (Precision* buffer) = 0;

    // Recover one heat-flux vector per element node into ELEMENT_NODAL storage
    virtual void compute_heat_flux(
        Field&       heat_flux,
        const Field& temperature
    ) = 0;
};

} // namespace fem::model
