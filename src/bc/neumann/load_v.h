/**
 * @file load_v.h
 * @brief Defines distributed body-force-density loads on structural elements.
 *
 * `VLoad` prescribes a force vector per physical volume. Structural elements
 * integrate that vector consistently with their interpolation and geometric
 * volume measure. Because the stored quantity is already a force density, no
 * material-density scaling is applied by this load type.
 *
 * @see Neumann
 * @see model::StructuralElement
 * @see cos::CoordinateSystem
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include "neumann.h"
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
struct VLoad : Neumann {
    // Types
    using Ptr = std::shared_ptr<VLoad>;

    // Nominal force-density components. NaN marks an omitted component.
    Vec3 values_ = {NAN, NAN, NAN};

    // Target structural elements over whose physical volume the field is
    // integrated.
    SPtr<model::ElementRegion> region_ = nullptr;

    // Construction
    VLoad() = default;
    ~VLoad() override = default;

    // Integrate N^T b over the selected element volumes without material-density
    // scaling and accumulate the equivalent nodal forces in the supplied RHS.
    void apply(model::ModelData& model_data, model::Field& rhs,
               Precision time, bool ignore_amplitude = false) override;

    // Diagnostics
    std::string str() const override;
};

} // namespace fem::bc
