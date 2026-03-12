/**
 * @file load_d.h
 * @brief Declares the distributed surface traction load.
 *
 * @see src/bc/load_d.cpp
 * @see src/bc/load.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "load.h"

#include "../data/region.h"

namespace fem {
namespace bc {

/**
 * @struct DLoad
 * @brief Distributed surface load with vector traction.
 *
 * Applies translational forces to the nodes of each surface patch referenced by
 * the associated region. If an orientation system is supplied, the traction
 * components are defined in that local frame and mapped to global coordinates
 * for every surface node.
 */
struct DLoad : public Load {
    using Ptr = std::shared_ptr<DLoad>; ///< Shared pointer alias for distributed loads.

    Vec3 values{NAN, NAN, NAN};                  ///< Surface traction components.
    SPtr<model::SurfaceRegion> region = nullptr; ///< Target surface region.

    DLoad() = default;
    ~DLoad() override = default;

    void apply(model::ModelData& model_data, model::Field& bc, Precision time) override;
    std::string str() const override;
};

} // namespace bc
} // namespace fem
