/**
 * @file load_c.h
 * @brief Declares concentrated nodal forces and moments.
 *
 * `CLoad` assigns up to three force and three moment components to every node
 * in a node region. Individual components may remain unspecified through
 * `NaN`, an optional coordinate system can define a local component basis, and
 * an optional amplitude scales the complete generalized load at analysis time.
 *
 * @see CLoad
 * @see Neumann
 * @see load_c.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "neumann.h"

#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Applies concentrated generalized Neumann loads to a set of nodes.
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

    // Assemble the concentrated generalized load into the nodal field.
    void apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude = false) override;

    std::string str() const override;
};

} // namespace fem::bc
