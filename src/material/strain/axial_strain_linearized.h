/**
 * @file axial_strain_linearized.h
 * @brief Declares the linearized axial strain measure.
 *
 * This type identifies a scalar axial strain as a small-strain quantity. Truss
 * elements use it for geometrically linear material evaluation with axial
 * Cauchy stress.
 *
 * @see axial_strain_linearized.cpp
 */

#pragma once

#include "axial_strain.h"

namespace fem {

// Small axial strain used with axial Cauchy stress
struct AxialStrainLinearized : AxialStrain {
    // Constructs a zero linearized axial strain
    AxialStrainLinearized() = default;

    // Constructs the strain directly from epsilon_xx
    explicit AxialStrainLinearized(Precision value);
};

} // namespace fem
