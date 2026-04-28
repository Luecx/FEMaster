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

    Vec3                       values_ = {NAN, NAN, NAN}; ///< Surface traction components.
    SPtr<model::SurfaceRegion> region_ = nullptr;          ///< Target surface region.

    /**
     * @brief Default constructor for delayed field assignment by collectors.
     */
    DLoad() = default;

    /**
     * @brief Defaulted virtual destructor.
     */
    ~DLoad() override = default;

    /**
     * @brief Integrates the traction over every target surface.
     *
     * @param model_data Model data used for positions and surface lookup.
     * @param bc Boundary-condition field receiving nodal forces.
     * @param time Current analysis time for amplitude scaling.
     */
    void apply(model::ModelData& model_data, model::Field& bc, Precision time) override;

    /**
     * @brief Returns a compact description of the distributed load.
     *
     * @return std::string Target region, values, orientation and amplitude.
     */
    std::string str() const override;
};
} // namespace bc
} // namespace fem
