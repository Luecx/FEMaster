/**
 * @file heat_flux.h
 * @brief Declares prescribed surface heat-flux boundary conditions.
 *
 * `HeatFlux` is a classical thermal Neumann condition. It prescribes heat flux
 * per unit boundary area and contributes only to the thermal right-hand side.
 * The existing three-component surface-vector integration path is reused; the
 * scalar thermal contribution is stored exclusively in component zero.
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
 * \f$q\f$, the consistent nodal source contribution is
 * \f$\mathbf{f}_q = \int_{\Gamma_q} \mathbf{N}^T q\,\mathrm{d}\Gamma\f$.
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

    /**
     * @brief Integrates prescribed heat flux into component zero of thermal RHS.
     *
     * @param model_data Compiled surface topology and reference geometry.
     * @param rhs Three-component nodal thermal load field. Only component zero
     *            is modified; components one and two remain untouched.
     * @param time Analysis time used for amplitude evaluation.
     * @param ignore_amplitude If true, omit amplitude scaling.
     * @param system_dof_ids Unused; prescribed heat flux does not modify LHS.
     * @param lhs Unused; prescribed heat flux does not modify LHS.
     */
    void apply(model::ModelData&       model_data,
               model::Field&           rhs,
               Precision               time,
               bool                    ignore_amplitude = false,
               const SystemDofIds*      system_dof_ids = nullptr,
               TripletList*             lhs = nullptr) override;

    std::string str() const override;
};

} // namespace fem::bc
