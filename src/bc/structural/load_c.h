/**
 * @file load_c.h
 * @brief Defines concentrated nodal forces and moments.
 *
 * `CLoad` is the direct nodal member of the structural load family. It stores up to
 * three force and three moment components and adds them to every node of a
 * target region after optional amplitude scaling and coordinate transformation.
 *
 * Unlike distributed loads, no finite-element quadrature is required because
 * the physical load is already defined at discrete nodal degrees of freedom.
 *
 * @see Condition
 * @see cos::CoordinateSystem
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include "../condition.h"
#include "../../constraints/types/equation.h"
#include "../../core/types_eig.h"
#include "../../data/field.h"
#include "../../cos/coordinate_system.h"
#include "../../data/region.h"

#include <cmath>
#include <memory>
#include <string>

namespace fem::bc {

/**
 * @brief Applies concentrated generalized loads to a node region.
 *
 * The six nominal components follow
 *
 *     [Fx, Fy, Fz, Mx, My, Mz].
 *
 * `NaN` marks an omitted component and therefore contributes zero. With scalar
 * amplitude `a(t)` and local-to-global basis `A(x)`, a prescribed local force
 * vector is assembled as
 *
 *     f_i <- f_i + a(t) A(x_i) f_local,
 *
 * while the moment triplet is transformed by the same basis and accumulated in
 * the rotational generalized DOFs. Without an orientation, `A = I`.
 */
struct CLoad : Condition {
    // Types
    using Ptr = std::shared_ptr<CLoad>;

    // Nominal generalized components. NaN omits a component.
    Vec6 values_ = {NAN, NAN, NAN, NAN, NAN, NAN};

    // Target nodes receiving the same nominal generalized load
    model::NodeRegion::Ptr region_ = nullptr;

    // Optional local basis in which vector components are prescribed
    cos::CoordinateSystem::Ptr orientation_ = nullptr;

    // Construction
    CLoad() = default;
    ~CLoad() override = default;

    // Assemble the generalized nodal loads defined by this condition into
    // the global load field. Each target node receives the prescribed force
    // and moment components [Fx, Fy, Fz, Mx, My, Mz].
    //
    // Initial and target step states scale this nominal load independently.
    //
    // Components marked as NaN are omitted. If a local coordinate system is
    // assigned, forces and moments are transformed into global coordinates
    // at each node before being accumulated in load. Dead loads evaluate the
    // basis at the reference position; follower loads evaluate it at the
    // current position and add the derivative of the rotating basis to
    // load_tangent. Constraints and the system matrix remain unchanged.

    void apply(
        model::ModelData&        model_data,
        model::Field&            load,
        TripletList*             load_tangent,
        constraint::Equations&   constraints,
        const model::SystemDofs& system_dof_ids,
        TripletList&             system_matrix,
        Precision                scale
    ) override;

    // Diagnostics
    std::string str() const override;
};

} // namespace fem::bc
