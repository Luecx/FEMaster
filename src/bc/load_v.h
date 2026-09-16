/**
 * @file load_v.h
 * @brief Declares distributed body-force loading on element regions.
 *
 * `VLoad` describes a force-density vector field that is integrated over the
 * volume, area or line measure represented by each selected structural element.
 * The element formulation multiplies the supplied vector by its shape functions
 * and the appropriate geometric measure before scattering equivalent nodal
 * forces. Material density is deliberately not applied. Optional orientation
 * and amplitude objects modify the nominal vector without changing the
 * element-specific integration logic.
 *
 * @see VLoad
 * @see Load
 * @see load_v.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "load.h"

#include "../data/region.h"

namespace fem {
namespace bc {

/**
 * @brief Applies a body-force density vector to structural elements.
 *
 * `NaN` components are treated as omitted and converted to zero. In the global
 * case the resulting vector field is spatially constant. With an orientation,
 * the nominal local vector is rotated at each integration point so coordinate
 * systems with position-dependent axes are supported.
 */
struct VLoad : public Load {
    // Shared ownership type for concrete volumetric-load references.
    using Ptr = std::shared_ptr<VLoad>;

    // Nominal body-force density vector, integrated directly over the element
    // measure without additional material-density scaling.
    Vec3 values_ = {NAN, NAN, NAN};

    // Element region whose structural elements receive the distributed load.
    SPtr<model::ElementRegion> region_ = nullptr;

    // Construct an empty body load for subsequent collector assignment.
    VLoad() = default;

    // Enable polymorphic destruction through `Load`.
    ~VLoad() override = default;

    // Build the global vector field, optionally rotate it at each integration
    // point and delegate consistent integration to each valid structural
    // element in `region_` without applying material density.
    void apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude = false) override;

    // Return the target element set, nominal components and any orientation or
    // amplitude in a compact diagnostic representation.
    std::string str() const override;
};
} // namespace bc
} // namespace fem
