/**
 * @file axial_strain_linearized.cpp
 * @brief Implements construction of linearized axial strains.
 *
 * The implementation initializes the common scalar `AxialStrain` storage while
 * preserving the distinct linearized strain type used for material dispatch.
 *
 * @see axial_strain_linearized.h
 */

#include "axial_strain_linearized.h"

namespace fem {

// Preserve the small-strain type while using the common scalar storage
AxialStrainLinearized::AxialStrainLinearized(Precision value)
    : AxialStrain(value) {}

} // namespace fem
