/**
 * @file load_v.h
 * @brief Declares the volumetric element load.
 *
 * @see src/bc/load_v.cpp
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
 * @struct VLoad
 * @brief Volumetric load applied to structural elements.
 *
 * Applies a body-force vector to each structural element in the associated
 * element region. For oriented loads the vector is taken in local coordinates
 * and rotated at the element centroid before being scattered to the nodes.
 */
struct VLoad : public Load {
    using Ptr = std::shared_ptr<VLoad>; ///< Shared pointer alias for volumetric loads.

    Vec3 values{NAN, NAN, NAN};                  ///< Body-force components.
    SPtr<model::ElementRegion> region = nullptr; ///< Target element region.

    VLoad() = default;
    ~VLoad() override = default;

    void apply(model::ModelData& model_data, model::Field& bc, Precision time) override;
    std::string str() const override;
};

} // namespace bc
} // namespace fem
