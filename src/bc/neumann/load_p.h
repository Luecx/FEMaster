/**
 * @file load_p.h
 * @brief Declares scalar pressure loading on surface regions.
 *
 * `PLoad` applies a uniform scalar pressure normal to every selected surface.
 * Surface geometry supplies the position-dependent normal and performs the
 * consistent integration into nodal forces. An optional amplitude scales the
 * pressure magnitude at the current analysis time.
 *
 * Pressure is a pure right-hand-side contribution. The optional system DOF map
 * and LHS triplet list carried by the common load interface are unused.
 *
 * @see PLoad
 * @see Neumann
 * @see load_p.cpp
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "neumann.h"

#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Integrates a uniform pressure as a structural load-side condition.
 *
 * The stored pressure is multiplied by the optional amplitude and applied in
 * the negative direction of the surface normal returned at each integration
 * point. Consequently, the vector direction follows each surface's geometric
 * orientation while the scalar magnitude remains uniform.
 */
struct PLoad : Neumann {
    using Ptr = std::shared_ptr<PLoad>;

    // Nominal scalar pressure magnitude.
    Precision pressure_ = NAN;

    // Surface region receiving the pressure.
    SPtr<model::SurfaceRegion> region_ = nullptr;

    PLoad() = default;
    ~PLoad() override = default;

    // Integrate pressure into the six-component structural RHS. The optional
    // LHS assembly objects are unused.
    void apply(model::ModelData&       model_data,
               model::Field&           rhs,
               Precision               time,
               bool                    ignore_amplitude = false,
               const SystemDofIds*      system_dof_ids = nullptr,
               TripletList*             lhs = nullptr) override;

    std::string str() const override;
};

} // namespace fem::bc
