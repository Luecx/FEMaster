/**
 * @file neumann.h
 * @brief Declares pure right-hand-side boundary conditions.
 *
 * Neumann conditions prescribe flux-like quantities that enter the weak finite-
 * element equilibrium only through the external right-hand side. They do not
 * add coefficients to the global operator and do not constrain primary nodal
 * unknowns. Structural and thermal specializations share this mathematical
 * behavior but target different nodal field layouts.
 *
 * The hierarchy intentionally remains separate from Robin conditions. A Robin
 * condition contains a state-dependent boundary term and therefore requires
 * simultaneous access to the right-hand side and system-operator assembly.
 *
 * @see Load
 * @see Neumann
 * @see StructuralNeumann
 * @see ThermalNeumann
 * @see ../robin/robin.h
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "../load.h"

#include "../../core/types_num.h"
#include "../../data/field.h"

#include <memory>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Common abstraction for boundary conditions contributing only to the RHS.
 *
 * In a weak finite-element formulation, a Neumann condition prescribes a flux or
 * generalized traction whose consistent nodal contribution can be assembled
 * without knowledge of the current primary unknown vector. Concrete conditions
 * therefore receive model data, a nodal target field and the current analysis
 * time, but no matrix or equation-assembly object.
 *
 * `amplitude_` and `orientation_` are inherited from `Load`. Concrete conditions
 * decide whether those shared modifiers are meaningful for their physical
 * quantity. `ignore_amplitude` is used when callers need the nominal spatial
 * load basis independent of its temporal history.
 */
struct Neumann : Load {
    // Shared ownership type used by load collectors and model-side dispatch.
    using Ptr = std::shared_ptr<Neumann>;

    // Enable destruction through the Neumann interface while retaining concrete
    // load ownership in shared-pointer collectors.
    virtual ~Neumann() = default;

    // Assemble the complete consistent nodal RHS contribution of this boundary
    // condition. The supplied field layout is determined by the structural or
    // thermal specialization and is validated by the concrete implementation.
    virtual void apply(model::ModelData& model_data,
                       model::Field&     rhs,
                       Precision         time,
                       bool              ignore_amplitude = false) = 0;
};

/**
 * @brief Marker base for mechanical Neumann conditions.
 *
 * Structural Neumann conditions assemble equivalent generalized nodal forces
 * into a NODE field with at least six components ordered as three translations
 * followed by three rotations. Concrete examples include concentrated forces,
 * pressure, distributed traction, body force and structural thermal expansion.
 *
 * The marker allows model-side collectors to distinguish mechanical load
 * assembly from scalar thermal heat-flow assembly without adding another virtual
 * dispatch method or duplicating common `Load` state.
 */
struct StructuralNeumann : Neumann {
    // Shared ownership type used when a caller needs to restrict a collection to
    // mechanical RHS-only boundary conditions.
    using Ptr = std::shared_ptr<StructuralNeumann>;

    ~StructuralNeumann() override = default;
};

/**
 * @brief Marker base for scalar thermal Neumann conditions.
 *
 * Thermal Neumann conditions prescribe heat flow into the stationary or
 * transient thermal balance. They assemble consistently into a one-component
 * NODE field whose sole component is the scalar thermal right-hand side. A
 * prescribed surface heat flux is the canonical implementation.
 *
 * Unlike structural thermal expansion loads, this hierarchy acts on the thermal
 * conduction system itself and never writes generalized mechanical forces.
 */
struct ThermalNeumann : Neumann {
    // Shared ownership type used by thermal load collectors and model-side
    // thermal boundary-condition dispatch.
    using Ptr = std::shared_ptr<ThermalNeumann>;

    ~ThermalNeumann() override = default;
};

} // namespace fem::bc
