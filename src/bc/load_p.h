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

    Precision                  pressure_ = NAN; ///< Magnitude of the applied pressure.
    SPtr<model::SurfaceRegion> region_   = nullptr; ///< Target surface region.

    /**
     * @brief Default constructor for delayed field assignment by collectors.
     */
    PLoad() = default;

    /**
     * @brief Defaulted virtual destructor.
     */
    ~PLoad() override = default;

    /**
     * @brief Integrates pressure over every target surface.
     *
     * @param model_data Model data used for positions and surface lookup.
     * @param bc Boundary-condition field receiving nodal forces.
     * @param time Current analysis time for amplitude scaling.
     */
    void apply(model::ModelData& model_data, model::Field& bc, Precision time) override;

    /**
     * @brief Returns a compact description of the pressure load.
     *
     * @return std::string Target region, pressure and amplitude.
     */
    std::string str() const override;
};
} // namespace bc
} // namespace fem
