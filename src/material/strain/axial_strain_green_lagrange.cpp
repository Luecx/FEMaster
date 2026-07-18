/**
 * @file axial_strain_green_lagrange.cpp
 * @brief Implements Green-Lagrange axial strain construction from stretch.
 *
 * The implementation converts an axial stretch ratio into the corresponding
 * Green-Lagrange strain and forwards direct scalar construction to the shared
 * `AxialStrain` storage.
 *
 * @see axial_strain_green_lagrange.h
 */

#include "axial_strain_green_lagrange.h"

namespace fem {

// Preserve the finite-strain type while using the common scalar storage
AxialStrainGreenLagrange::AxialStrainGreenLagrange(Precision value)
    : AxialStrain(value) {}

// Convert axial stretch into E_xx = 0.5 * (lambda^2 - 1)
AxialStrainGreenLagrange AxialStrainGreenLagrange::from_stretch(Precision stretch) {
    return AxialStrainGreenLagrange(
        Precision(0.5) * (stretch * stretch - Precision(1))
    );
}

} // namespace fem
