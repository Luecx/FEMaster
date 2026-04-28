/**
 * @file load_c.h
 * @brief Declares the concentrated nodal load boundary condition.
 *
 * @see src/bc/load_c.cpp
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
 * @struct CLoad
 * @brief Concentrated nodal load.
 *
 * Applies up to six generalized force/moment components to each node in the
 * associated node region. When an orientation coordinate system is provided,
 * the components are interpreted in the local basis and rotated into the
 * global frame prior to assembly.
 */
struct CLoad : public Load {
    using Ptr = std::shared_ptr<CLoad>; ///< Shared pointer alias for concentrated loads.

    Vec6                    values_ = {NAN, NAN, NAN, NAN, NAN, NAN}; ///< Generalized load vector (Fx,Fy,Fz,Mx,My,Mz).
    SPtr<model::NodeRegion> region_ = nullptr;                         ///< Target node region.

    /**
     * @brief Default constructor for delayed field assignment by collectors.
     */
    CLoad() = default;

    /**
     * @brief Defaulted virtual destructor.
     */
    ~CLoad() override = default;

    /**
     * @brief Scatters concentrated forces and moments to each target node.
     *
     * @param model_data Model data used for node positions and local orientation lookup.
     * @param bc Boundary-condition field receiving force and moment components.
     * @param time Current analysis time for amplitude scaling.
     */
    void apply(model::ModelData& model_data, model::Field& bc, Precision time) override;

    /**
     * @brief Returns a compact description of the concentrated load.
     *
     * @return std::string Target region, values, orientation and amplitude.
     */
    std::string str() const override;
};
} // namespace bc
} // namespace fem
