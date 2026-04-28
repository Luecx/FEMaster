/**
 * @file load_t.h
 * @brief Declares the thermal load boundary condition.
 *
 * @see src/bc/load_t.cpp
 * @see src/bc/load.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "load.h"

namespace fem {
namespace bc {
/**
 * @struct TLoad
 * @brief Thermal load derived from a temperature field.
 *
 * Generates equivalent nodal forces by feeding a temperature field to each
 * structural element.
 */
struct TLoad : public Load {
    using Ptr = std::shared_ptr<TLoad>; ///< Shared pointer alias for thermal loads.

    SPtr<model::Field> temp_field_ = nullptr; ///< Temperature field reference.
    Precision          ref_temp_   = NAN;     ///< Reference temperature for zero load.

    /**
     * @brief Default constructor for delayed field assignment by collectors.
     */
    TLoad() = default;

    /**
     * @brief Defaulted virtual destructor.
     */
    ~TLoad() override = default;

    /**
     * @brief Applies thermal expansion loads to all structural elements.
     *
     * @param model_data Model data that provides the structural elements.
     * @param bc Boundary-condition field receiving equivalent nodal forces.
     * @param time Current analysis time. Thermal loads currently ignore it.
     */
    void apply(model::ModelData& model_data, model::Field& bc, Precision time) override;

    /**
     * @brief Returns a compact description of the thermal load.
     *
     * @return std::string Temperature field name and reference temperature.
     */
    std::string str() const override;
};
} // namespace bc
} // namespace fem
