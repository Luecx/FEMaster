/**
 * @file elasticity.cpp
 * @brief Implements default behavior of the constitutive interface.
 *
 * The base implementation declares every kinematic formulation unsupported and
 * provides history integration only for stateless laws. Stateful constitutive
 * models must override `integrate()` so a source state is never modified in
 * place and the complete updated history is written to a separate target row.
 *
 * @see Elasticity
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#include "elasticity.h"

#include "strain/axial_strain_green_lagrange.h"
#include "strain/axial_strain_linearized.h"
#include "strain/beam_generalized_strain.h"
#include "strain/shell_material_strain_green_lagrange.h"
#include "strain/shell_material_strain_linearized.h"
#include "strain/volume_strain_green_lagrange.h"
#include "strain/volume_strain_linearized.h"
#include "stress/axial_stress_cauchy.h"
#include "stress/axial_stress_pk2.h"
#include "stress/beam_stress_resultants.h"
#include "stress/shell_material_stress_cauchy.h"
#include "stress/shell_material_stress_pk2.h"
#include "stress/volume_stress_cauchy.h"
#include "stress/volume_stress_pk2.h"

namespace fem::material {

bool Elasticity::supports_axial_linearized() const { return false; }
bool Elasticity::supports_axial_green_lagrange() const { return false; }
bool Elasticity::supports_volume_linearized() const { return false; }
bool Elasticity::supports_volume_green_lagrange() const { return false; }
bool Elasticity::supports_beam_resultants() const { return false; }
bool Elasticity::supports_shell_integration_linearized() const { return false; }
bool Elasticity::supports_shell_integration_green_lagrange() const { return false; }

Index Elasticity::state_size() const { return 0; }

void Elasticity::initialize_state(Precision* state) const {
    (void) state;
}

void Elasticity::evaluate(const AxialStrainLinearized& strain,
                          const Precision*             state,
                          AxialStressCauchy&           stress,
                          Precision&                   tangent) const {
    (void) strain; (void) state; (void) stress; (void) tangent;
    logging::error(false, "Elasticity model does not support linearized axial evaluation");
}

void Elasticity::evaluate(const AxialStrainGreenLagrange& strain,
                          const Precision*                  state,
                          AxialStressPK2&                 stress,
                          Precision&                     tangent) const {
    (void) strain; (void) state; (void) stress; (void) tangent;
    logging::error(false, "Elasticity model does not support Green-Lagrange axial evaluation");
}

void Elasticity::evaluate(const VolumeStrainLinearized& strain,
                          const Precision*              state,
                          VolumeStressCauchy&           stress,
                          Mat6&                         tangent) const {
    (void) strain; (void) state; (void) stress; (void) tangent;
    logging::error(false, "Elasticity model does not support linearized volume evaluation");
}

void Elasticity::evaluate(const VolumeStrainGreenLagrange& strain,
                          const Precision*                   state,
                          VolumeStressPK2&                 stress,
                          Mat6&                           tangent) const {
    (void) strain; (void) state; (void) stress; (void) tangent;
    logging::error(false, "Elasticity model does not support Green-Lagrange volume evaluation");
}

void Elasticity::evaluate(const BeamGeneralizedStrain& strain,
                          const Precision*             state,
                          BeamStressResultants&        resultants,
                          Mat6&                        tangent) const {
    (void) strain; (void) state; (void) resultants; (void) tangent;
    logging::error(false, "Elasticity model does not support beam-resultant evaluation");
}

void Elasticity::evaluate(const ShellMaterialStrainLinearized& strain,
                          const Precision*                     state,
                          ShellMaterialStressCauchy&            stress,
                          Mat5&                               tangent) const {
    (void) strain; (void) state; (void) stress; (void) tangent;
    logging::error(false, "Elasticity model does not support linearized shell evaluation");
}

void Elasticity::evaluate(const ShellMaterialStrainGreenLagrange& strain,
                          const Precision*                          state,
                          ShellMaterialStressPK2&                 stress,
                          Mat5&                                  tangent) const {
    (void) strain; (void) state; (void) stress; (void) tangent;
    logging::error(false, "Elasticity model does not support Green-Lagrange shell evaluation");
}

namespace {

void require_stateless_integration(const Elasticity& elasticity,
                                   const Precision* state,
                                   Precision* target_state) {
    logging::error(elasticity.state_size() == 0,
        "Stateful elasticity model must override constitutive integration");
    (void) state;
    (void) target_state;
}

} // namespace

void Elasticity::integrate(const AxialStrainLinearized& strain,
                           const Precision*             state,
                           Precision*                   target_state,
                           AxialStressCauchy&           stress,
                           Precision&                   tangent) const {
    require_stateless_integration(*this, state, target_state);
    evaluate(strain, state, stress, tangent);
}

void Elasticity::integrate(const AxialStrainGreenLagrange& strain,
                           const Precision*                  state,
                           Precision*                      target_state,
                           AxialStressPK2&                 stress,
                           Precision&                     tangent) const {
    require_stateless_integration(*this, state, target_state);
    evaluate(strain, state, stress, tangent);
}

void Elasticity::integrate(const VolumeStrainLinearized& strain,
                           const Precision*              state,
                           Precision*                    target_state,
                           VolumeStressCauchy&           stress,
                           Mat6&                         tangent) const {
    require_stateless_integration(*this, state, target_state);
    evaluate(strain, state, stress, tangent);
}

void Elasticity::integrate(const VolumeStrainGreenLagrange& strain,
                           const Precision*                   state,
                           Precision*                       target_state,
                           VolumeStressPK2&                 stress,
                           Mat6&                           tangent) const {
    require_stateless_integration(*this, state, target_state);
    evaluate(strain, state, stress, tangent);
}

void Elasticity::integrate(const BeamGeneralizedStrain& strain,
                           const Precision*             state,
                           Precision*                   target_state,
                           BeamStressResultants&        resultants,
                           Mat6&                        tangent) const {
    require_stateless_integration(*this, state, target_state);
    evaluate(strain, state, resultants, tangent);
}

void Elasticity::integrate(const ShellMaterialStrainLinearized& strain,
                           const Precision*                     state,
                           Precision*                           target_state,
                           ShellMaterialStressCauchy&            stress,
                           Mat5&                               tangent) const {
    require_stateless_integration(*this, state, target_state);
    evaluate(strain, state, stress, tangent);
}

void Elasticity::integrate(const ShellMaterialStrainGreenLagrange& strain,
                           const Precision*                          state,
                           Precision*                              target_state,
                           ShellMaterialStressPK2&                 stress,
                           Mat5&                                  tangent) const {
    require_stateless_integration(*this, state, target_state);
    evaluate(strain, state, stress, tangent);
}

} // namespace fem::material
