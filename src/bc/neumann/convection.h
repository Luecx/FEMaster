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
 * \f$q = h(T_\infty - T)\f$.
 * After finite-element discretization and rearrangement, convection contributes
 *
 * \f$\mathbf{K}_h = \int_\Gamma h\,\mathbf{N}^T\mathbf{N}\,\mathrm{d}\Gamma\f$
 *
 * to the left-hand side and
 *
 * \f$\mathbf{f}_h = \int_\Gamma hT_\infty\,\mathbf{N}^T\,\mathrm{d}\Gamma\f$
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

    /**
     * @brief Assembles both convection RHS and boundary-matrix contributions.
     *
     * The RHS contribution is integrated through the existing three-component
     * surface-vector integration path and stored only in component zero. The LHS
     * contribution is assembled directly into sparse triplets using the supplied
     * thermal DOF map.
     *
     * @param model_data Compiled surface topology and reference geometry.
     * @param rhs Three-component nodal thermal load field. Only component zero
     *            receives the ambient-temperature source term.
     * @param time Analysis time used for amplitude evaluation.
     * @param ignore_amplitude If true, omit amplitude scaling.
     * @param system_dof_ids Thermal node-to-system DOF map. Must be non-null.
     * @param lhs Sparse triplet list receiving the convection boundary matrix.
     *            Must be non-null.
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
