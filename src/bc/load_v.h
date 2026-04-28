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

    Vec3                       values_ = {NAN, NAN, NAN}; ///< Body-force components.
    SPtr<model::ElementRegion> region_ = nullptr;          ///< Target element region.

    /**
     * @brief Default constructor for delayed field assignment by collectors.
     */
    VLoad() = default;

    /**
     * @brief Defaulted virtual destructor.
     */
    ~VLoad() override = default;

    /**
     * @brief Integrates the body-force vector over every structural element.
     *
     * @param model_data Model data used for element and position lookup.
     * @param bc Boundary-condition field receiving nodal forces.
     * @param time Current analysis time for amplitude scaling.
     */
    void apply(model::ModelData& model_data, model::Field& bc, Precision time) override;

    /**
     * @brief Returns a compact description of the volumetric load.
     *
     * @return std::string Target region, values, orientation and amplitude.
     */
    std::string str() const override;
};
} // namespace bc
} // namespace fem
