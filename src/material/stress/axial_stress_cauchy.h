#pragma once

#include "axial_stress.h"

namespace fem {

struct AxialStressCauchy : AxialStress {
    AxialStressCauchy() = default;
    explicit AxialStressCauchy(Precision value);
};

} // namespace fem
