/**
 * @file load_v.h
 * @brief Declares distributed body-force loading over element regions.
 *
 * `VLoad` defines a vector body force per unit mass. Structural elements
 * integrate the field with material density and distribute the resulting force
 * consistently to their nodes. Optional orientation and amplitude modifiers are
 * inherited from `Neumann`.
 *
 * Body force contributes only to the structural right-hand side. The optional
 * system DOF map and LHS triplet list from the common interface are unused.
 *
 * @see VLoad
 * @see Neumann
 * @see load_v.cpp
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "neumann.h"

#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Integrates density-scaled body forces over selected elements.
 */
struct VLoad : Neumann {
    using Ptr = std::shared_ptr<VLoad>;

    // Nominal body-force components per unit mass.
    Vec3 values_ = {NAN, NAN, NAN};

    // Element region receiving the distributed load.
    SPtr<model::ElementRegion> region_ = nullptr;

    VLoad() = default;
    ~VLoad() override = default;

    /**
     * @brief Integrates the density-scaled body-force field into nodal RHS.
     *
     * @param model_data Compiled element topology, materials and positions.
     * @param rhs Structural nodal right-hand-side field.
     * @param time Analysis time used for amplitude evaluation.
     * @param ignore_amplitude If true, omit amplitude scaling.
     * @param system_dof_ids Unused; body force does not modify LHS.
     * @param lhs Unused; body force does not modify LHS.
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
