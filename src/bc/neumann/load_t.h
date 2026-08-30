/**
 * @file load_t.h
 * @brief Declares structural thermal-expansion loading from a nodal temperature field.
 *
 * `TLoad` does not prescribe equivalent nodal forces directly. It validates a
 * scalar nodal temperature field and delegates thermal-strain evaluation to
 * every structural element together with the stress-free reference temperature.
 *
 * The resulting equivalent thermal force is a pure structural right-hand-side
 * contribution. The optional system DOF map and LHS triplet list exposed by the
 * common load interface are therefore unused by this condition.
 *
 * @see TLoad
 * @see Neumann
 * @see load_t.cpp
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "neumann.h"

namespace fem::bc {

/**
 * @brief Applies element-specific thermal expansion loads to the model.
 *
 * The condition owns no spatial target region because the supplied nodal
 * temperature field already covers the compiled model. Every structural
 * element extracts its local temperatures, evaluates thermal strain relative
 * to `ref_temp_` and scatters the corresponding equivalent nodal force into the
 * structural RHS. This is a mechanical thermal-expansion load, not a heat-flow
 * boundary condition.
 */
struct TLoad : Neumann {
    using Ptr = std::shared_ptr<TLoad>;

    // Scalar nodal temperature field consumed by structural elements.
    SPtr<model::Field> temp_field_ = nullptr;

    // Stress-free reference temperature.
    Precision ref_temp_ = NAN;

    TLoad() = default;
    ~TLoad() override = default;

    // Convert the prescribed temperature field into the six-component
    // structural RHS. Time, amplitude and the optional LHS objects are unused.
    void apply(model::ModelData&       model_data,
               model::Field&           rhs,
               Precision               time,
               bool                    ignore_amplitude = false,
               const SystemDofIds*      system_dof_ids = nullptr,
               TripletList*             lhs = nullptr) override;

    std::string str() const override;
};

} // namespace fem::bc
