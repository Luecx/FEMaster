/**
 * @file axial_strain.cpp
 * @brief Implements axial Green-Lagrange strain storage and construction.
 *
 * This material kinematics implementation stores the scalar reference strain
 * and constructs E_xx = 0.5 * (lambda^2 - 1) from axial stretch. Component access
 * also supports increments and variations used by linearized element response.
 * Constitutive evaluation and material history remain caller responsibilities.
 *
 * @see axial_strain.h
 *
 * @author Finn Eggers
 * @date 09.10.2026
 */

#include "axial_strain.h"

namespace fem {

// Initialize the common scalar storage
AxialStrain::AxialStrain(Precision value)
    : value_(value) {}

// Access the only axial component through the component interface
Precision& AxialStrain::operator[](Component) {
    return value_;
}

Precision AxialStrain::operator[](Component) const {
    return value_;
}

// Access the scalar directly for constitutive evaluation
Precision AxialStrain::value() const {
    return value_;
}

Precision& AxialStrain::value() {
    return value_;
}

/**
 * @brief Constructs axial Green-Lagrange strain from the stretch ratio.
 *
 * The stretch lambda is the current length divided by the reference length.
 * This constructs the total strain E_xx = 0.5 * (lambda^2 - 1), rather than
 * a linearized increment about a supplied base state. No material state changes.
 *
 * @param stretch Axial stretch ratio supplied by the element kinematics.
 * @return Total axial Green-Lagrange strain in the reference direction.
 */
AxialStrain AxialStrain::from_stretch(Precision stretch) {
    // Construct the axial component of E = 0.5 * (F^T F - I).
    return AxialStrain(Precision(0.5) * (stretch * stretch - Precision(1)));
}

} // namespace fem
