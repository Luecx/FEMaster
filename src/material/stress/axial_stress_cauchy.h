/**
 * @file axial_stress_cauchy.h
 * @brief Declares the axial Cauchy stress measure.
 *
 * The scalar represents true axial stress in the current configuration. It is
 * used with linearized axial strain and by material paths that return spatial
 * stress directly.
 *
 * @see axial_stress_cauchy.cpp
 */

#pragma once

#include "axial_stress.h"

namespace fem {

// Axial true stress defined in the current configuration
struct AxialStressCauchy : AxialStress {
    // Constructs a zero axial Cauchy stress
    AxialStressCauchy() = default;

    // Constructs directly from the axial Cauchy stress value
    explicit AxialStressCauchy(Precision value);
};

} // namespace fem
