/**
 * @file convection.h
 * @brief Defines linear thermal convection as a Mixed boundary condition.
 *
 * Convection couples the unknown surface temperature to a prescribed ambient
 * temperature through a film coefficient. The weak form contributes both an
 * ambient heat-flow source and a symmetric boundary operator, so the condition
 * belongs to the Mixed algebraic category.
 *
 * @see Mixed
 * @see ThermalCondition
 * @see model::SurfaceInterface
 *
 * @author Finn Eggers
 * @date 18.09.2026
 */

#pragma once

#include "mixed.h"
#include "../thermal.h"
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
struct Convection : Mixed, ThermalCondition {
    // Types
    using Ptr = std::shared_ptr<Convection>;

    // Compiled surface region exchanging heat with the ambient environment
    model::SurfaceRegion::Ptr region_ = nullptr;

    // Nominal film coefficient and prescribed ambient temperature
    Precision film_coefficient_    = Precision(0);
    Precision ambient_temperature_ = Precision(0);

    // Assemble q_h = integral h T_inf N^T dGamma into the scalar thermal RHS
    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) override;

    // Assemble K_h = integral h N N^T dGamma into active thermal system triplets
    void apply_matrix(model::ModelData&   model_data,
                      const SystemDofIds& system_dof_ids,
                      TripletList&        matrix,
                      Precision           time,
                      bool                ignore_amplitude = false) override;

    // Diagnostics
    std::string str() const override;

private:
    // Evaluate and validate the amplitude-scaled film coefficient used by both
    // RHS and operator assembly.
    Precision effective_film_coefficient(Precision time, bool ignore_amplitude) const;
};

} // namespace fem::bc
