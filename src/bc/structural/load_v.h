/**
 * @file load_v.h
 * @brief Defines distributed body-force-density loads on structural elements.
 *
 * `VLoad` prescribes a force vector per physical volume. Structural elements
 * integrate that vector consistently with their interpolation and geometric
 * volume measure. Because the stored quantity is already a force density, no
 * material-density scaling is applied by this load type.
 *
 * @see Condition
 * @see model::StructuralElement
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
 * @brief Applies a distributed force density to an element region.
 *
 * For one structural element, the consistent nodal contribution is
 *
 *     f_e = integral_Omega_e N^T b dOmega,
 *
 * where `b` has units of force per volume. If the stored components are defined
 * in a local basis, the physical vector field becomes
 *
 *     b(x) = a(t) A(x) b_local.
 *
 * `NaN` entries denote omitted vector components and contribute zero. The
 * element integration is explicitly performed without density scaling; inertia
 * loads use a separate formulation where mass density is part of the measure.
 */
struct VLoad : Condition {
    // Types
    using Ptr = std::shared_ptr<VLoad>;

    // Target and initial body-force density components [bx, by, bz].
    // NaN in values_ omits a component, while NaN in values_start_ indicates
    // an unspecified initial value, interpreted as zero during interpolation.
    // The effective force density is assembled per unit physical volume.
    Vec3 values_       = {NAN, NAN, NAN};
    Vec3 values_start_ = {NAN, NAN, NAN};

    // Target structural elements over whose physical volumes the load is integrated
    model::ElementRegion::Ptr region_ = nullptr;

    // Optional local basis in which vector components are prescribed
    cos::CoordinateSystem::Ptr orientation_ = nullptr;

    // Construction
    VLoad() = default;
    ~VLoad() override = default;

    // Integrate the prescribed body-force density over every structural element
    // in region_ and accumulate the consistent equivalent forces in the RHS.
    //
    // Without an amplitude, values_start_ and values_ are interpolated using
    // step_progress. With an amplitude, the target is scaled by its value at
    // time instead. ignore_amplitude applies the nominal target directly.
    //
    // NaN target components are omitted. An optional local coordinate system
    // is transformed to global coordinates at each integration point.
    // Integration uses the physical volume measure without material-density
    // scaling because this load is already a force per unit volume.
    //
    // Only rhs is modified; equations, system_dof_ids and matrix are unchanged.
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
