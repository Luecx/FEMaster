/**
 * @file load.h
 * @brief Defines the common right-hand-side assembly capability of load-like conditions.
 *
 * `Load` is not a boundary-condition category. It is the shared assembly
 * interface used by `Neumann` and `Mixed` conditions that contribute to the
 * global right-hand side. Concrete conditions retain their mathematical
 * classification through the boundary-condition hierarchy while reusing the
 * common orientation and amplitude state defined here.
 *
 * For a finite-element equilibrium problem written schematically as
 *
 *     R(u) = f_int(u) - f_ext = 0,
 *
 * `Load::apply()` contributes terms to the external generalized load field
 * `f_ext`. Mixed conditions may additionally contribute an operator term through
 * their own matrix-assembly interface.
 *
 * @see Neumann
 * @see Mixed
 * @see LoadCollector
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include "amplitude.h"
#include "../core/printable.h"
#include "../core/types_cls.h"
#include "../core/types_eig.h"
#include "../cos/coordinate_system.h"
#include "../data/field.h"

#include <memory>
#include <string>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Shared right-hand-side assembly interface for load-like conditions.
 *
 * A load maps its physical definition to equivalent generalized nodal forces.
 * The exact mapping is owned by the concrete condition: concentrated loads add
 * nodal values directly, surface and volume loads perform consistent finite-
 * element integration, and thermal or inertial loads construct equivalent
 * structural forces from element kinematics and material data.
 *
 * `orientation_` defines the basis in which vector-valued load components are
 * prescribed. Concrete conditions evaluate that basis at the geometrically
 * appropriate position before transforming the components to global
 * coordinates. `amplitude_` optionally supplies the scalar time history
 * multiplying the nominal load definition.
 *
 * The interface intentionally does not inherit from `BoundaryCondition` so it
 * can be shared by the independent `Neumann` and `Mixed` categories without
 * introducing an artificial inheritance relation between them.
 */
struct Load : Printable {
    // Types
    using Ptr = std::shared_ptr<Load>;

    // Optional load modifiers shared by RHS-capable conditions
    cos::CoordinateSystem::Ptr orientation_ = nullptr;
    Amplitude::Ptr             amplitude_   = nullptr;

    // Polymorphic destruction
    virtual ~Load() = default;

    // Assemble the equivalent generalized nodal load into `rhs`. The concrete
    // condition owns interpolation, geometric integration and coordinate
    // transformation. `ignore_amplitude` evaluates the nominal load with unit
    // temporal scale.
    virtual void apply(model::ModelData& model_data,
                       model::Field&     rhs,
                       Precision         time,
                       bool              ignore_amplitude = false) = 0;

    // Diagnostics
    std::string str() const override = 0;
};

} // namespace fem::bc
