/**
 * @file beam_generalized_strain.cpp
 * @brief Implements storage and component access for generalized beam strains.
 *
 * This file maps named beam strain components onto the underlying six-entry
 * vector and provides direct vector access for element and section assembly.
 *
 * @see beam_generalized_strain.h
 */

#include "beam_generalized_strain.h"

namespace fem {

// Initialize all generalized beam components in formulation order
BeamGeneralizedStrain::BeamGeneralizedStrain(const Vec6& values)
    : values_(values) {}

// Map named components onto the underlying vector
Precision& BeamGeneralizedStrain::operator[](Component component) {
    return values_(static_cast<int>(component));
}

Precision BeamGeneralizedStrain::operator[](Component component) const {
    return values_(static_cast<int>(component));
}

// Expose the complete state for section evaluation and assembly
const Vec6& BeamGeneralizedStrain::values() const {
    return values_;
}

Vec6& BeamGeneralizedStrain::values() {
    return values_;
}

} // namespace fem
