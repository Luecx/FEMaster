#pragma once

#include "axial_stress.h"

namespace fem {

struct AxialStressPK2 : AxialStress {
    AxialStressPK2() = default;
    explicit AxialStressPK2(Precision value);
};

} // namespace fem
