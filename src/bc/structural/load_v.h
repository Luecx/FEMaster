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

    // Nominal generalized components. NaN omits a component.
    Vec3 values_ = {NAN, NAN, NAN};

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
    // Each step evaluates this nominal force density using its initial or
    // target amplitude factor, without modifying the physical definition.
    //
    // NaN target components are omitted. An optional local coordinate system
    // is transformed to global coordinates at each integration point.
    // Integration uses the physical volume measure without material-density
    // scaling because this load is already a force per unit volume.
    //
    // Only rhs is modified; equations, system_dof_ids and matrix are unchanged.
    // Construct an independent physical condition for a new step.
    Condition::Ptr clone() const override { return std::make_shared<VLoad>(*this); }

    void apply(
        model::ModelData&      model_data,
        model::Field&          rhs,
        constraint::Equations& equations,
        const SystemDofIds&    system_dof_ids,
        TripletList&           matrix,
        Precision              time,
        bool                   ignore_amplitude = false,
        Precision              step_progress = Precision(1),
        bool initial = false,
        Precision previous_end_time = Precision(0),
        Precision end_time = Precision(1)
    ) override;

    // Diagnostics
    std::string str() const override;
};

} // namespace fem::bc
