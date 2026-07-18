/**
 * @file axial_strain_green_lagrange.h
 * @brief Declares the Green-Lagrange axial strain measure.
 *
 * The stored scalar represents `E = 0.5 * (lambda^2 - 1)`. It is used by
 * geometrically nonlinear truss elements when evaluating finite-strain
 * materials together with an energetically conjugate PK2 stress.
 *
 * @see axial_strain_green_lagrange.cpp
 */

#pragma once

#include "axial_strain.h"

namespace fem {

// Finite axial strain energetically conjugate to axial PK2 stress
struct AxialStrainGreenLagrange : AxialStrain {
    // Constructs a zero Green-Lagrange axial strain
    AxialStrainGreenLagrange() = default;

    // Constructs the strain directly from E_xx
    explicit AxialStrainGreenLagrange(Precision value);

    // Computes E_xx = 0.5 * (stretch^2 - 1) from the axial stretch
    static AxialStrainGreenLagrange from_stretch(Precision stretch);
};

} // namespace fem
