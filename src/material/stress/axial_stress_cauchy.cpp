/**
 * @file axial_stress_cauchy.cpp
 * @brief Implements construction of axial Cauchy stress.
 *
 * The implementation reuses the scalar `AxialStress` storage while preserving
 * the spatial Cauchy-stress type used for material dispatch.
 *
 * @see axial_stress_cauchy.h
 */

#include "axial_stress_cauchy.h"

namespace fem {

// Preserve the Cauchy-stress type while using common scalar storage
AxialStressCauchy::AxialStressCauchy(Precision value)
    : AxialStress(value) {}

} // namespace fem
