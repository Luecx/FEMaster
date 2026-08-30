/**
 * @file convection.h
 * @brief Declares linear thermal convection boundary conditions.
 *
 * `Convection` is mathematically a Robin boundary condition, but it participates
 * in FEMaster's common load-side `Neumann` hierarchy so all load-like boundary
 * conditions use one assembly interface. In one `apply()` call it contributes
 * both the ambient-temperature source term to the thermal right-hand side and
 * the temperature-dependent film term to the global left-hand-side operator.
 *
 * @see Convection
 * @see Neumann
 * @see convection.cpp
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "neumann.h"

namespace fem::bc {

/**
 * @brief Applies a linear convection law to a surface region.
 *
 * The boundary law is
 *
 *     q = h (T_inf - T).
 *
 * After finite-element discretization and rearrangement, convection contributes
 *
 *     K_h = integral_Gamma h N^T N dGamma
 *
 * to the left-hand side and
 *
 *     f_h = integral_Gamma h T_inf N^T dGamma
 *
 * to the right-hand side. Both contributions are assembled by the single common
 * `apply()` method. Consequently `system_dof_ids` and `lhs` are mandatory when
 * this condition is evaluated.
 */
struct Convection : Neumann {
    using Ptr = std::shared_ptr<Convection>;

    // Boundary surfaces exposed to the surrounding medium.
    model::SurfaceRegion::Ptr region_ = nullptr;

    // Film coefficient h relating temperature difference to outward heat flow.
    Precision film_coefficient_ = Precision(0);

    // Prescribed ambient/sink temperature T_inf.
    Precision ambient_temperature_ = Precision(0);

    Convection() = default;
    ~Convection() override = default;

    // Assemble the scalar ambient-temperature source and the matching boundary
    // matrix using the supplied thermal system mapping.
    void apply(model::ModelData&       model_data,
               model::Field&           rhs,
               Precision               time,
               bool                    ignore_amplitude = false,
               const SystemDofIds*      system_dof_ids = nullptr,
               TripletList*             lhs = nullptr) override;

    std::string str() const override;
};

} // namespace fem::bc
