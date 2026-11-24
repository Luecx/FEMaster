/**
 * @file lamina_elasticity.h
 * @brief Declares orthotropic plane-stress elasticity for shell laminae.
 *
 * A lamina is defined by the six constants `E1`, `E2`, `nu12`, `G12`, `G13`
 * and `G23` in its local material basis, matching the Abaqus `TYPE=LAMINA`
 * definition. The material provides only the shell response because the
 * through-thickness normal stiffness is not part of the definition.
 *
 * @see Elasticity
 * @see CompositeShellSection
 *
 * @author Finn Eggers
 * @date 24.11.2026
 */

#pragma once

#include "elasticity.h"

namespace fem::material {

/**
 * @brief Stateless orthotropic plane-stress elasticity for shell sections.
 *
 * The in-plane response uses the orthotropic plane-stress tangent with the
 * reciprocal ratio `nu21 = nu12 E2/E1`. Transverse shear uses the independent
 * moduli `G13` and `G23`. Solid, beam and truss evaluations are rejected by the
 * `Elasticity` base class.
 */
struct LaminaElasticity : Elasticity {
    // Young's moduli along the in-plane principal material directions
    Precision E1;
    Precision E2;

    // Major in-plane Poisson ratio: strain in 2 for loading in 1
    Precision nu12;

    // Engineering shear moduli in the principal material planes
    Precision G12;
    Precision G13;
    Precision G23;

    // Construct the lamina and validate positive moduli and a positive-definite
    // plane-stress tangent.
    LaminaElasticity(Precision E1,
                     Precision E2,
                     Precision nu12,
                     Precision G12,
                     Precision G13,
                     Precision G23);

    // Finite-strain shell response returning PK2 components work-conjugate to
    // the five supplied Green-Lagrange strain components. State remains unchanged.
    void evaluate(const ShellMaterialStrain&              strain,
                  const Precision*                        old_state,
                  Precision*                              new_state,
                  ShellMaterialStressPK2&                 stress,
                  Mat5*                                   tangent = nullptr) const override;
};

} // namespace fem::material
