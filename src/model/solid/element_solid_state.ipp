/**
 * @file element_solid_state.ipp
 * @brief Implements solid constitutive source-to-target state integration.
 *
 * Read-only stress recovery remains in `evaluate_material()`. The routines in
 * this file are used only when a nonlinear equilibrium evaluation is allowed to
 * advance constitutive history from the accepted source field into the working
 * target field.
 *
 * @author Finn Eggers
 * @date 31.08.2026
 */

#pragma once

#include "../../section/section_solid.h"

namespace fem::model {

/**
 * Integrates a linearized solid material state from source to target history.
 *
 * The section performs the material-basis transformation and constitutive
 * integration. Optional element stiffness scaling is applied afterwards to the
 * returned stress and tangent, exactly as in the read-only evaluation path.
 */
template<Index N>
void SolidElement<N>::integrate_material(
    Precision                     r,
    Precision                     s,
    Precision                     t,
    const VolumeStrainLinearized& global_strain,
    const Precision*              state,
    Precision*                    target_state,
    VolumeStressCauchy&           global_stress,
    Mat6&                         global_tangent
) {
    get_section()->integrate(
        material_position_reference(r, s, t),
        additional_material_rotation(),
        global_strain,
        state,
        target_state,
        global_stress,
        global_tangent
    );

    const Precision scaling = element_stiffness_scale();
    global_stress.voigt() *= scaling;
    global_tangent        *= scaling;
}

/**
 * Integrates a Total-Lagrangian solid material state from source to target.
 *
 * Green-Lagrange strain and PK2 stress remain in the global reference basis.
 * The target row is fully written by the material integration and may therefore
 * be promoted atomically by the nonlinear state manager after convergence.
 */
template<Index N>
void SolidElement<N>::integrate_material(
    Precision                        r,
    Precision                        s,
    Precision                        t,
    const VolumeStrainGreenLagrange& global_strain,
    const Precision*                 state,
    Precision*                       target_state,
    VolumeStressPK2&                 global_stress,
    Mat6&                            global_tangent
) {
    get_section()->integrate(
        material_position_reference(r, s, t),
        additional_material_rotation(),
        global_strain,
        state,
        target_state,
        global_stress,
        global_tangent
    );

    const Precision scaling = element_stiffness_scale();
    global_stress.voigt() *= scaling;
    global_tangent        *= scaling;
}

} // namespace fem::model
