/**
 * @file load_t.h
 * @brief Declares structural thermal-expansion loading.
 */

#pragma once

#include "neumann.h"

namespace fem::bc {

struct TLoad : StructuralNeumann {
    using Ptr = std::shared_ptr<TLoad>;

    SPtr<model::Field> temp_field_ = nullptr;
    Precision ref_temp_ = NAN;

    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) override;

    std::string str() const override;
};

} // namespace fem::bc
