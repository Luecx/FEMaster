/**
 * @file beam_stress_resultants.cpp
 * @brief Implements storage access for generalized beam stress resultants.
 *
 * The six-entry vector is passed between beam constitutive evaluation and
 * element assembly without changing the beam formulation's component order.
 *
 * @see beam_stress_resultants.h
 */

#include "beam_stress_resultants.h"

namespace fem {

// Initialize all generalized beam forces and moments
BeamStressResultants::BeamStressResultants(const Vec6& values)
    : values_(values) {}

// Expose the complete resultant state for beam assembly
const Vec6& BeamStressResultants::values() const {
    return values_;
}

Vec6& BeamStressResultants::values() {
    return values_;
}

} // namespace fem
