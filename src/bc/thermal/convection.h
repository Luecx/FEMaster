/**
 * @file convection.h
 * @brief Defines linear thermal convection as a thermal boundary condition.
 *
 * Convection couples the unknown surface temperature to a prescribed ambient
 * temperature through a film coefficient. The weak form contributes both an
 * ambient heat-flow source and a symmetric boundary operator. Those two
 * operations are exposed directly by this concrete thermal condition without an
 * intermediate algebraic base class.
 *
 * @see Condition
 * @see model::SurfaceInterface
 *
 * @author Finn Eggers
 * @date 18.09.2026
 */

#pragma once

#include "../condition.h"
#include "../../constraints/types/equation.h"
#include "../../core/types_eig.h"
#include "../../data/field.h"
#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Applies Newton cooling on a compiled surface region.
 *
 * With outward conductive flux written as
 *
 *     q_n = h (T - T_inf),
 *
 * the finite-element weak form contributes
 *
 *     K_h = integral_Gamma h N N^T dGamma
 *
 * and
 *
 *     q_h = integral_Gamma h T_inf N^T dGamma.
 *
 * The resulting stationary thermal balance contains `K_T + K_h` on the left
 * and the ambient source `q_h` on the right. The optional amplitude scales the
 * film coefficient `h`, not the ambient temperature.
 */
struct Convection : Condition {
    // Types
    using Ptr = std::shared_ptr<Convection>;

    // Compiled surfaces exchanging heat with the surrounding environment
    model::SurfaceRegion::Ptr region_ = nullptr;

    // Film coefficient h and ambient temperature T_inf. An assigned amplitude
    // scales h, and therefore both the ambient source and boundary operator.
    Precision film_coefficient_    = Precision(0);
    Precision ambient_temperature_ = Precision(0);

    // Assemble the complete Newton cooling contribution into the thermal system.
    //
    // The temperature-dependent flux q_n = h * (T - T_inf) contributes
    //
    //     K_h = integral_Gamma h N N^T dGamma
    //     q_h = integral_Gamma h T_inf N^T dGamma
    //
    // The ambient source q_h is accumulated in rhs, while K_h is assembled
    // into matrix using system_dof_ids. The surface integrals use the reference
    // geometry and consistent finite-element shape functions.
    //
    // An assigned amplitude scales the film coefficient h at physical time;
    // ignore_amplitude uses the nominal coefficient. There is currently no
    // start/target history for convection, so step_progress is unused.
    // Constraint equations remain unchanged.
    void apply(
        model::ModelData&      model_data,
        model::Field&          rhs,
        constraint::Equations& equations,
        const SystemDofIds&    system_dof_ids,
        TripletList&           matrix,
        Precision              time,
        bool                   ignore_amplitude = false,
        Precision              step_progress   = Precision(1)
    ) override;

    // Assemble only the ambient heat source integral q_h into the scalar RHS.
    // This method does not modify the film boundary operator.
    void apply_rhs(
        model::ModelData& model_data,
        model::Field&     rhs,
        Precision         time,
        bool              ignore_amplitude = false
    );

    // Assemble only K_h into the sparse thermal system triplets. The nodal
    // temperature DOFs are mapped to equation indices via system_dof_ids.
    void apply_matrix(
        model::ModelData&   model_data,
        const SystemDofIds& system_dof_ids,
        TripletList&        matrix,
        Precision           time,
        bool                ignore_amplitude = false
    );

    // Diagnostics
    std::string str() const override;

private:
    // Evaluate and validate the amplitude-scaled film coefficient used by both
    // RHS and operator assembly.
    Precision effective_film_coefficient(Precision time, bool ignore_amplitude) const;
};

} // namespace fem::bc
