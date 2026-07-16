#include "elasticity.h"

#include "strain/axial_strain_green_lagrange.h"
#include "strain/axial_strain_linearized.h"
#include "strain/beam_generalized_strain.h"
#include "strain/shell_generalized_strain.h"
#include "strain/shell_material_strain_green_lagrange.h"
#include "strain/shell_material_strain_linearized.h"
#include "strain/volume_strain_green_lagrange.h"
#include "strain/volume_strain_linearized.h"
#include "stress/axial_stress_cauchy.h"
#include "stress/axial_stress_pk2.h"
#include "stress/beam_stress_resultants.h"
#include "stress/shell_material_stress_cauchy.h"
#include "stress/shell_material_stress_pk2.h"
#include "stress/shell_stress_resultants.h"
#include "stress/volume_stress_cauchy.h"
#include "stress/volume_stress_pk2.h"

namespace fem::material {

bool Elasticity::supports_axial_linearized() const {
    return false;
}

bool Elasticity::supports_axial_green_lagrange() const {
    return false;
}

bool Elasticity::supports_volume_linearized() const {
    return false;
}

bool Elasticity::supports_volume_green_lagrange() const {
    return false;
}

bool Elasticity::supports_beam_resultants() const {
    return false;
}

bool Elasticity::supports_shell_resultants() const {
    return false;
}

bool Elasticity::supports_shell_integration_linearized() const {
    return false;
}

bool Elasticity::supports_shell_integration_green_lagrange() const {
    return false;
}

Index Elasticity::state_size() const {
    return 0;
}

void Elasticity::initialize_state(Precision* state) const {
    (void) state;
}

void Elasticity::evaluate(const AxialStrainLinearized& strain,
                          const Precision*             state_old,
                          Precision*                   state_new,
                          AxialStressCauchy&           stress,
                          Precision&                   tangent) const {
    (void) strain;
    (void) state_old;
    (void) state_new;
    (void) stress;
    (void) tangent;
    logging::error(false, "Elasticity model does not support linearized axial evaluation");
}

void Elasticity::evaluate(const AxialStrainGreenLagrange& strain,
                          const Precision*                state_old,
                          Precision*                      state_new,
                          AxialStressPK2&                 stress,
                          Precision&                      tangent) const {
    (void) strain;
    (void) state_old;
    (void) state_new;
    (void) stress;
    (void) tangent;
    logging::error(false, "Elasticity model does not support Green-Lagrange axial evaluation");
}

void Elasticity::evaluate(const VolumeStrainLinearized& strain,
                          const Precision*              state_old,
                          Precision*                    state_new,
                          VolumeStressCauchy&           stress,
                          Mat6&                         tangent) const {
    (void) strain;
    (void) state_old;
    (void) state_new;
    (void) stress;
    (void) tangent;
    logging::error(false, "Elasticity model does not support linearized volume evaluation");
}

void Elasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                          const Precision*                 state_old,
                          Precision*                       state_new,
                          VolumeStressPK2&                 stress,
                          Mat6&                            tangent) const {
    (void) strain;
    (void) state_old;
    (void) state_new;
    (void) stress;
    (void) tangent;
    logging::error(false, "Elasticity model does not support Green-Lagrange volume evaluation");
}

void Elasticity::evaluate(const BeamGeneralizedStrain& strain,
                          const Precision*             state_old,
                          Precision*                   state_new,
                          BeamStressResultants&        resultants,
                          Mat6&                        tangent) const {
    (void) strain;
    (void) state_old;
    (void) state_new;
    (void) resultants;
    (void) tangent;
    logging::error(false, "Elasticity model does not support beam-resultant evaluation");
}

void Elasticity::evaluate(const ShellGeneralizedStrain& strain,
                          Precision                     thickness,
                          const Precision*              state_old,
                          Precision*                    state_new,
                          ShellStressResultants&        resultants,
                          Mat8&                         tangent) const {
    (void) strain;
    (void) thickness;
    (void) state_old;
    (void) state_new;
    (void) resultants;
    (void) tangent;
    logging::error(false, "Elasticity model does not support shell-resultant evaluation");
}

void Elasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                          const Precision*                     state_old,
                          Precision*                           state_new,
                          ShellMaterialStressCauchy&            stress,
                          Mat5&                                 tangent) const {
    (void) strain;
    (void) state_old;
    (void) state_new;
    (void) stress;
    (void) tangent;
    logging::error(false, "Elasticity model does not support linearized shell integration");
}

void Elasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                          const Precision*                        state_old,
                          Precision*                              state_new,
                          ShellMaterialStressPK2&                 stress,
                          Mat5&                                   tangent) const {
    (void) strain;
    (void) state_old;
    (void) state_new;
    (void) stress;
    (void) tangent;
    logging::error(false, "Elasticity model does not support Green-Lagrange shell integration");
}

} // namespace fem::material
