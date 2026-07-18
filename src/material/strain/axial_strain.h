/**
 * @file axial_strain.h
 * @brief Declares the base type for scalar axial strain measures.
 *
 * `AxialStrain` stores the single normal-strain component required by truss
 * material models. Derived types distinguish linearized and Green-Lagrange
 * kinematics while sharing the same scalar storage and component interface.
 *
 * @see axial_strain.cpp
 */

#pragma once

#include "../../core/types_num.h"

namespace fem {

// Common storage interface for a scalar axial strain
struct AxialStrain {
    // Named access to the single available strain component
    enum class Component : int {
        XX = 0
    };

    // Constructs a zero axial strain
    AxialStrain() = default;

    // Constructs an axial strain from its scalar value
    explicit AxialStrain(Precision value);

    // Returns the selected component by mutable or constant access
    Precision& operator[](Component component);
    Precision  operator[](Component component) const;

    // Returns the scalar strain value by constant or mutable access
    [[nodiscard]] Precision  value() const;
    [[nodiscard]] Precision& value();

protected:
    // Stored axial strain component
    Precision value_{};
};

} // namespace fem
