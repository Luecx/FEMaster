/**
 * @file lamina_elasticity.cpp
 * @brief Implements orthotropic plane-stress elasticity for shell laminae.
 *
 * The shell material tangent is ordered `[11,22,12,13,23]`. Its in-plane block
 * is the orthotropic plane-stress reduction and the transverse shear entries are
 * the prescribed moduli `G13` and `G23`.
 *
 * @see LaminaElasticity
 *
 * @author Finn Eggers
 * @date 24.11.2026
 */

#include "lamina_elasticity.h"

#include "../core/logging.h"
#include "strain/shell_material_strain.h"
#include "stress/shell_material_stress_cauchy.h"
#include "stress/shell_material_stress_pk2.h"

namespace fem::material {

/**
 * Constructs lamina elasticity from in-plane engineering constants.
 *
 * Stability of the plane-stress tangent requires positive moduli and
 * `nu12^2 E2/E1 < 1`.
 *
 * @param E1 Young's modulus along material direction 1.
 * @param E2 Young's modulus along material direction 2.
 * @param nu12 Poisson ratio for loading in direction 1 and contraction in 2.
 * @param G12 Engineering shear modulus in the 1-2 plane.
 * @param G13 Engineering shear modulus in the 1-3 plane.
 * @param G23 Engineering shear modulus in the 2-3 plane.
 */
LaminaElasticity::LaminaElasticity(Precision E1,
                                   Precision E2,
                                   Precision nu12,
                                   Precision G12,
                                   Precision G13,
                                   Precision G23)
    : E1  (E1),
      E2  (E2),
      nu12(nu12),
      G12 (G12),
      G13 (G13),
      G23 (G23) {
    logging::error(E1 > Precision(0) && E2 > Precision(0),
        "LaminaElasticity: E1 and E2 must be positive");
    logging::error(G12 > Precision(0) && G13 > Precision(0) && G23 > Precision(0),
        "LaminaElasticity: G12, G13 and G23 must be positive");
    logging::error(nu12 * nu12 * E2 / E1 < Precision(1),
        "LaminaElasticity: nu12^2 E2/E1 must be smaller than 1");
}

/**
 * Evaluates lamina shell PK2 stress from Green-Lagrange strain.
 *
 * @param strain Five-component Green-Lagrange material strain.
 * @param old_state Unused input material-point state row.
 * @param new_state Unused output material-point state row.
 * @param stress Shell second Piola-Kirchhoff stress.
 * @param tangent Optional reduced material derivative.
 */
void LaminaElasticity::evaluate(const ShellMaterialStrain&              strain,
                                const Precision*                        old_state,
                                Precision*                              new_state,
                                ShellMaterialStressPK2&                 stress,
                                Mat5*                                   tangent) const {
    (void) old_state;
    (void) new_state;

    const Precision nu21  = nu12 * E2 / E1;
    const Precision denom = Precision(1) - nu12 * nu21;

    Mat5 material_tangent = Mat5::Zero();
    material_tangent(0, 0) = E1 / denom;
    material_tangent(0, 1) = nu12 * E2 / denom;
    material_tangent(1, 0) = nu21 * E1 / denom;
    material_tangent(1, 1) = E2 / denom;
    material_tangent(2, 2) = G12;
    material_tangent(3, 3) = G13;
    material_tangent(4, 4) = G23;

    stress.values() = material_tangent * strain.values();

    if (tangent != nullptr) {
        *tangent = material_tangent;
    }
}

} // namespace fem::material
