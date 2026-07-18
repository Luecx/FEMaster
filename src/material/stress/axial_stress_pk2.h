/**
 * @file axial_stress_pk2.h
 * @brief Declares the axial second Piola-Kirchhoff stress measure.
 *
 * The scalar is defined in the reference configuration and is energetically
 * conjugate to Green-Lagrange axial strain. Geometrically nonlinear truss
 * material models use this type for their constitutive response.
 *
 * @see axial_stress_pk2.cpp
 */

#pragma once

#include "axial_stress.h"

namespace fem {

// Axial PK2 stress defined in the reference configuration
struct AxialStressPK2 : AxialStress {
    // Constructs a zero axial PK2 stress
    AxialStressPK2() = default;

    // Constructs directly from the axial PK2 stress value
    explicit AxialStressPK2(Precision value);
};

} // namespace fem
