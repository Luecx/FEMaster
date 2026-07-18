/**
 * @file axial_stress.cpp
 * @brief Implements storage and value access for scalar axial stresses.
 *
 * The shared implementation is independent of the actual stress measure;
 * Cauchy and PK2 semantics are carried by their derived types.
 *
 * @see axial_stress.h
 */

#include "axial_stress.h"

namespace fem {

// Initialize the common scalar storage
AxialStress::AxialStress(Precision value)
    : value_(value) {}

// Access the scalar directly for constitutive evaluation
Precision AxialStress::value() const {
    return value_;
}

Precision& AxialStress::value() {
    return value_;
}

} // namespace fem
