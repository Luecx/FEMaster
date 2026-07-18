/**
 * @file axial_stress.h
 * @brief Declares the base type for scalar axial stress measures.
 *
 * `AxialStress` stores the single normal-stress component returned by truss
 * material models. Derived types distinguish Cauchy and PK2 stress while
 * sharing the same scalar storage interface.
 *
 * @see axial_stress.cpp
 */

#pragma once

#include "../../core/types_num.h"

namespace fem {

// Common storage interface for a scalar axial stress
struct AxialStress {
    // Constructs a zero axial stress
    AxialStress() = default;

    // Constructs an axial stress from its scalar value
    explicit AxialStress(Precision value);

    // Returns the scalar stress value by constant or mutable access
    [[nodiscard]] Precision  value() const;
    [[nodiscard]] Precision& value();

protected:
    // Stored axial stress component
    Precision value_{};
};

} // namespace fem
