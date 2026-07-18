/**
 * @file beam_generalized_strain.h
 * @brief Declares generalized axial, shear and curvature strains for beams.
 *
 * The six-component state contains axial strain, two transverse shear strains
 * and three curvatures. Beam elements use it as the section-level kinematic
 * state from which axial force, shear forces, torsion and bending moments are
 * evaluated.
 *
 * @see beam_generalized_strain.cpp
 */

#pragma once

#include "../../core/types_eig.h"

namespace fem {

// Section-level kinematic state of a beam
struct BeamGeneralizedStrain {
    // Component order used by the internal six-entry vector
    enum class Component : Index {
        EpsilonX = 0,
        GammaY   = 1,
        GammaZ   = 2,
        KappaX   = 3,
        KappaY   = 4,
        KappaZ   = 5
    };

    // Constructs a zero generalized beam strain
    BeamGeneralizedStrain() = default;

    // Constructs the state in the order defined by Component
    explicit BeamGeneralizedStrain(const Vec6& values);

    // Returns a named component by mutable or constant access
    Precision& operator[](Component component);
    Precision  operator[](Component component) const;

    // Returns all generalized components by constant or mutable access
    [[nodiscard]] const Vec6& values() const;
    [[nodiscard]] Vec6&       values();

private:
    // [epsilon_x, gamma_y, gamma_z, kappa_x, kappa_y, kappa_z]
    Vec6 values_{Vec6::Zero()};
};

} // namespace fem
