/**
 * @file load_d.h
 * @brief Declares distributed vector tractions on surface regions.
 *
 * `DLoad` represents a traction vector per unit surface area. Each selected
 * surface performs shape-function integration and scatters the consistent nodal
 * force. The vector may use a local coordinate system and an amplitude.
 *
 * The condition contributes only to the structural right-hand side. The
 * optional system DOF map and LHS triplet list from the common load interface
 * are accepted for uniformity and intentionally ignored.
 *
 * @see DLoad
 * @see Neumann
 * @see load_d.cpp
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "neumann.h"

#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Integrates a prescribed traction vector over selected surfaces.
 *
 * `values_` represents force per unit physical area. For each surface, the
 * condition evaluates the possibly position-dependent coordinate-system basis,
 * applies amplitude scaling and uses the surface shape functions to construct a
 * consistent generalized nodal load. Rotational RHS components remain zero.
 */
struct DLoad : Neumann {
    using Ptr = std::shared_ptr<DLoad>;

    // Nominal surface-traction components. NaN denotes an omitted component.
    Vec3 values_ = {NAN, NAN, NAN};

    // Surface region receiving the traction.
    SPtr<model::SurfaceRegion> region_ = nullptr;

    DLoad() = default;
    ~DLoad() override = default;

    // Integrate the traction into the six-component structural RHS. The
    // optional LHS assembly objects are unused.
    void apply(model::ModelData&       model_data,
               model::Field&           rhs,
               Precision               time,
               bool                    ignore_amplitude = false,
               const SystemDofIds*      system_dof_ids = nullptr,
               TripletList*             lhs = nullptr) override;

    std::string str() const override;
};

} // namespace fem::bc
