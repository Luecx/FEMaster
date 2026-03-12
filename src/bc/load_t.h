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

    SPtr<model::Field> temp_field = nullptr; ///< Temperature field reference.
    Precision ref_temp{NAN};                 ///< Reference temperature for zero load.

    TLoad() = default;
    ~TLoad() override = default;

    void apply(model::ModelData& model_data, model::Field& bc, Precision time) override;
    std::string str() const override;
};

} // namespace bc
} // namespace fem
