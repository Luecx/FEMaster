/**
 * @file axial_strain.h
 * @brief Declares axial Green-Lagrange strain and its increments.
 *
 * The material kinematics subsystem uses AxialStrain for the scalar reference
 * strain E_xx = 0.5 * (lambda^2 - 1), work-conjugate to axial PK2 stress.
 * Elements construct total strains or their linearized increments; constitutive
 * evaluation and ownership of material history remain with the material model
 * and its caller. This file provides storage, component access and construction
 * from the axial stretch ratio.
 *
 * @see axial_strain.cpp
 *
 * @author Finn Eggers
 * @date 09.10.2026
 */

#pragma once

#include "../../core/types_num.h"

namespace fem {

/**
 * @brief Scalar axial Green-Lagrange strain in the reference member direction.
 *
 * Stores E_xx, or an increment/variation dE_xx, with identical component access.
 * The caller distinguishes total strains from increments and supplies a total
 * state to nonlinear constitutive evaluation; increments enter through dS/dE.
 * Linearization about the undeformed state gives the infinitesimal axial strain.
 * The scalar is owned by value, initialized to zero, and carries no material
 * history. AxialStressPK2 is its work-conjugate stress representation.
 */
struct AxialStrain {
    // Named access to the single available strain component
    enum class Component : Index {
        XX = 0
    };

    // Constructs a zero axial strain
    AxialStrain() = default;

    // Constructs an axial strain from its scalar value
    explicit AxialStrain(Precision value);

    // Computes E_xx = 0.5 * (stretch^2 - 1) from the axial stretch
    static AxialStrain from_stretch(Precision stretch);

    // Returns the selected component by mutable or constant access
    Precision& operator[](Component component);
    Precision  operator[](Component component) const;

    // Returns the scalar strain value by constant or mutable access
    [[nodiscard]] Precision  value() const;
    [[nodiscard]] Precision& value();

private:
    // Total E_xx or its increment/variation
    Precision value_{};
};

} // namespace fem
