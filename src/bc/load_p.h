/**
 * @file load_p.h
 * @brief Declares the uniform surface pressure load.
 *
 * @see src/bc/load_p.cpp
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
 * @struct PLoad
 * @brief Uniform surface pressure load.
 *
 * Distributes a scalar pressure to the nodes of each surface in the associated
 * region.
 */
struct PLoad : public Load {
    using Ptr = std::shared_ptr<PLoad>; ///< Shared pointer alias for pressure loads.

    Precision pressure{NAN};               ///< Magnitude of the applied pressure.
    SPtr<model::SurfaceRegion> region = nullptr; ///< Target surface region.

    PLoad() = default;
    ~PLoad() override = default;

    void apply(model::ModelData& model_data, model::Field& bc, Precision time) override;
    std::string str() const override;
};

} // namespace bc
} // namespace fem
