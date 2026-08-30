/**
 * @file load_c.h
 * @brief Declares concentrated nodal forces and moments.
 *
 * `CLoad` assigns up to three force and three moment components to every node
 * in a node region. Individual components may remain unspecified through
 * `NaN`, an optional coordinate system can define a local component basis, and
 * an optional amplitude scales the complete generalized load at analysis time.
 *
 * `CLoad` is a pure right-hand-side condition. The common `Neumann::apply()`
 * interface also exposes optional system-matrix assembly objects for conditions
 * such as convection; `CLoad` deliberately ignores those objects.
 *
 * @see CLoad
 * @see Neumann
 * @see load_c.cpp
 * @author Finn Eggers
 * @date 30.08.2026
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
 * each node and rotated into the global frame before assembly.
 */
struct CLoad : Neumann {
    using Ptr = std::shared_ptr<CLoad>;

    // Nominal generalized nodal load [Fx, Fy, Fz, Mx, My, Mz].
    Vec6 values_ = {NAN, NAN, NAN, NAN, NAN, NAN};

    // Node region receiving the same generalized load.
    SPtr<model::NodeRegion> region_ = nullptr;

    CLoad() = default;
    ~CLoad() override = default;

    // Add the concentrated generalized load to the six-component structural
    // RHS. The optional LHS assembly objects are unused.
    void apply(model::ModelData&       model_data,
               model::Field&           rhs,
               Precision               time,
               bool                    ignore_amplitude = false,
               const SystemDofIds*      system_dof_ids = nullptr,
               TripletList*             lhs = nullptr) override;

    std::string str() const override;
};

} // namespace fem::bc
