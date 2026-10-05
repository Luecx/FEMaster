/**
 * @file elasticity.cpp
 * @brief Implements default behavior of the elastic constitutive interface.
 *
 * The base implementation defines a stateless material-point layout and rejects
 * constitutive overloads that a concrete elasticity model has not implemented.
 *
 * @see Elasticity
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#include "elasticity.h"

#include "strain/axial_strain_green_lagrange.h"
#include "strain/shell_material_strain_green_lagrange.h"
#include "strain/volume_strain_green_lagrange.h"
#include "stress/axial_stress_cauchy.h"
#include "stress/axial_stress_pk2.h"
#include "stress/shell_material_stress_cauchy.h"
#include "stress/shell_material_stress_pk2.h"
#include "stress/volume_stress_cauchy.h"
#include "stress/volume_stress_pk2.h"

namespace fem::material {

Index Elasticity::state_size() const {
    return 0;
}

void Elasticity::initialize_state(Precision* state) const {
    (void) state;
}

void Elasticity::evaluate(const AxialStrainGreenLagrange& strain,
                          const Precision*                old_state,
                          Precision*                      new_state,
                          AxialStressPK2&                 stress,
                          Precision*                      tangent) const {
    (void) strain;
    (void) old_state;
    (void) new_state;
    (void) stress;
    (void) tangent;

    logging::error(false,
        "Elasticity model does not support Green-Lagrange axial evaluation");
}

void Elasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                          const Precision*                 old_state,
                          Precision*                       new_state,
                          VolumeStressPK2&                 stress,
                          Mat6*                            tangent) const {
    (void) strain;
    (void) old_state;
    (void) new_state;
    (void) stress;
    (void) tangent;

    logging::error(false,
        "Elasticity model does not support Green-Lagrange volume evaluation");
}

void Elasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                          const Precision*                        old_state,
                          Precision*                              new_state,
                          ShellMaterialStressPK2&                 stress,
                          Mat5*                                   tangent) const {
    (void) strain;
    (void) old_state;
    (void) new_state;
    (void) stress;
    (void) tangent;

    logging::error(false,
        "Elasticity model does not support Green-Lagrange shell integration");
}

} // namespace fem::material
