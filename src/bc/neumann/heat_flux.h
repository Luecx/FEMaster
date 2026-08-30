/**
 * @file heat_flux.h
 * @brief Declares prescribed surface heat-flux boundary conditions.
 *
 * `HeatFlux` is a classical thermal Neumann condition. It prescribes heat flux
 * per unit boundary area and contributes only to a one-component thermal
 * right-hand side. Surface interpolation distributes the scalar boundary input
 * consistently to the connected temperature DOFs.
 *
 * @see HeatFlux
 * @see Neumann
 * @see heat_flux.cpp
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "neumann.h"

namespace fem::bc {

/**
 * @brief Applies a uniform prescribed heat flux to a surface region.
 *
 * Positive `heat_flux_` denotes heat entering the model. For a constant flux
 * q, the consistent nodal source contribution is
 *
 *     f_q = integral_Gamma_q N^T q dGamma.
 *
 * No temperature-dependent term is introduced, so the global left-hand-side
 * operator is not modified.
 */
struct HeatFlux : Neumann {
    using Ptr = std::shared_ptr<HeatFlux>;

    // Boundary surfaces on which the prescribed heat flux acts.
    model::SurfaceRegion::Ptr region_ = nullptr;

    // Nominal heat flux per unit area. Positive values enter the model.
    Precision heat_flux_ = Precision(0);

    HeatFlux() = default;
    ~HeatFlux() override = default;

    // Integrate the prescribed flux consistently into a scalar nodal thermal
    // right-hand side. Heat flux does not modify the system matrix.
    void apply(model::ModelData&       model_data,
               model::Field&           rhs,
               Precision               time,
               bool                    ignore_amplitude = false,
               const SystemDofIds*      system_dof_ids = nullptr,
               TripletList*             lhs = nullptr) override;

    std::string str() const override;
};

} // namespace fem::bc
