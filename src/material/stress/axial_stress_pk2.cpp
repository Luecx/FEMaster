/**
 * @file axial_stress_pk2.cpp
 * @brief Implements construction of axial PK2 stress.
 *
 * The implementation reuses the scalar `AxialStress` storage while preserving
 * the reference-configuration PK2 type used for material dispatch.
 *
 * @see axial_stress_pk2.h
 */

#include "axial_stress_pk2.h"

namespace fem {

// Preserve the PK2-stress type while using common scalar storage
AxialStressPK2::AxialStressPK2(Precision value)
    : AxialStress(value) {}

} // namespace fem
