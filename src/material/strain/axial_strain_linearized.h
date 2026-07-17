#pragma once

#include "axial_strain.h"

namespace fem {

struct AxialStrainLinearized : AxialStrain {
    AxialStrainLinearized() = default;
    explicit AxialStrainLinearized(Precision value);
};

} // namespace fem
