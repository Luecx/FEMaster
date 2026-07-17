#include "axial_strain_green_lagrange.h"

namespace fem {

AxialStrainGreenLagrange::AxialStrainGreenLagrange(Precision value)
    : AxialStrain(value) {}

AxialStrainGreenLagrange AxialStrainGreenLagrange::from_stretch(Precision stretch) {
    return AxialStrainGreenLagrange(
        Precision(0.5) * (stretch * stretch - Precision(1))
    );
}

} // namespace fem
