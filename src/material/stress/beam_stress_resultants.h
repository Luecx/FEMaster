/**
 * @file beam_stress_resultants.h
 * @brief Declares the generalized section-force state of a beam.
 *
 * The six entries store the force and moment resultants associated with the six
 * generalized beam strains. Beam constitutive paths return this state together
 * with a six-by-six section tangent.
 *
 * @see beam_stress_resultants.cpp
 */

#pragma once

#include "../../core/types_eig.h"

namespace fem {

// Generalized beam forces and moments in section coordinates
struct BeamStressResultants {
    // Constructs zero beam stress resultants
    BeamStressResultants() = default;

    // Constructs the resultants in the beam formulation's component order
    explicit BeamStressResultants(const Vec6& values);

    // Returns all section resultants by constant or mutable access
    [[nodiscard]] const Vec6& values() const;
    [[nodiscard]] Vec6&       values();

private:
    // Six generalized forces energetically conjugate to BeamGeneralizedStrain
    Vec6 values_{Vec6::Zero()};
};

} // namespace fem
