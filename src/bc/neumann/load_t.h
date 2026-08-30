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
 */
struct TLoad : Neumann {
    using Ptr = std::shared_ptr<TLoad>;

    // Scalar nodal temperature field consumed by structural elements.
    SPtr<model::Field> temp_field_ = nullptr;

    // Stress-free reference temperature.
    Precision ref_temp_ = NAN;

    TLoad() = default;
    ~TLoad() override = default;

    /**
     * @brief Converts the prescribed temperature field into equivalent nodal RHS.
     *
     * @param model_data Compiled structural element topology and section data.
     * @param rhs Structural nodal right-hand-side field receiving thermal force.
     * @param time Unused analysis time retained by the common interface.
     * @param ignore_amplitude Unused; `TLoad` currently has no amplitude scaling.
     * @param system_dof_ids Unused; structural thermal expansion does not modify LHS.
     * @param lhs Unused; structural thermal expansion does not modify LHS.
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
