#pragma once

#include "axial_strain.h"

namespace fem {

struct AxialStrainGreenLagrange : AxialStrain {
    AxialStrainGreenLagrange() = default;
    explicit AxialStrainGreenLagrange(Precision value);

    static AxialStrainGreenLagrange from_stretch(Precision stretch);
};

} // namespace fem
