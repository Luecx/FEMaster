/**
 * @file load_c.h
 * @brief Declares concentrated nodal forces and moments.
 *
 * `CLoad` assigns up to three force and three moment components to every node
 * in a node region. Individual components may remain unspecified through
 * `NaN`, an optional coordinate system can define a local component basis, and
 * an optional amplitude scales the complete generalized load at analysis time.
 * Assembly and diagnostic formatting are implemented in `load_c.cpp`.
 *
 * As a structural Neumann condition, the load contributes only to the global
 * six-component structural right-hand side and never modifies the system
 * operator or the structural constraint equations.
 *
 * @see CLoad
 * @see StructuralNeumann
 * @see Load
 * @see load_c.cpp
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "neumann.h"
#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Applies concentrated generalized loads to a set of nodes.
 *
 * The first three entries of `values_` represent translational forces and the
 * final three represent moments. `NaN` marks an omitted component. When an
 * orientation is assigned, both triplets are interpreted in its local basis at
 * each node and rotated into the global frame before they are accumulated in
 * the structural right-hand-side field.
 */
struct CLoad : StructuralNeumann {
    // Shared ownership type for code that needs the concrete concentrated-load
    // interface instead of the generic Neumann or Load pointer type.
    using Ptr = std::shared_ptr<CLoad>;

    // Nominal generalized nodal load ordered as [Fx, Fy, Fz, Mx, My, Mz]. NaN
    // denotes an unspecified component that contributes nothing during assembly.
    Vec6 values_ = {NAN, NAN, NAN, NAN, NAN, NAN};

    // Node region receiving the same generalized load at every contained node.
    SPtr<model::NodeRegion> region_ = nullptr;

    // Sanitize omitted components, apply the optional amplitude and assemble
    // forces and moments for every node in `region_`. Local components are
    // rotated with `orientation_` at the corresponding nodal position.
    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) override;

    // Describe the target node set, six nominal components and any assigned
    // orientation or amplitude in a compact diagnostic string.
    std::string str() const override;
};

} // namespace fem::bc
