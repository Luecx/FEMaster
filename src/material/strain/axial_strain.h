#pragma once

#include "../../core/types_num.h"

namespace fem {

struct AxialStrain {
    AxialStrain() = default;
    explicit AxialStrain(Precision value);

    [[nodiscard]] Precision  value() const;
    [[nodiscard]] Precision& value();

protected:
    Precision value_{};
};

} // namespace fem
