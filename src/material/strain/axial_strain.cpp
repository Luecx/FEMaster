/**
 * @file axial_strain.cpp
 * @brief Implements storage and component access for axial strains.
 *
 * This file provides the scalar constructor and the common component and value
 * access used by both linearized and finite-strain truss material evaluations.
 * It contains no constitutive behavior itself.
 *
 * @see axial_strain.h
 */

#include "axial_strain.h"

namespace fem {

// Initialize the common scalar storage
AxialStrain::AxialStrain(Precision value)
    : value_(value) {}

// Access the only axial component through the component interface
Precision& AxialStrain::operator[](Component component) {
    return value_;
}

Precision AxialStrain::operator[](Component component) const {
    return value_;
}

// Access the scalar directly for constitutive evaluation
Precision AxialStrain::value() const {
    return value_;
}

Precision& AxialStrain::value() {
    return value_;
}

} // namespace fem
