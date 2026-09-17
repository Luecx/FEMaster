/**
 * @file load_t.h
 * @brief Defines structural equivalent loading from a prescribed temperature field.
 *
 * `TLoad` is a structural load, not a thermal primary-variable boundary
 * condition. It references a scalar nodal temperature field and a stress-free
 * reference temperature, then delegates the conversion of thermal expansion to
 * equivalent nodal forces to each structural element formulation.
 *
 * This condition belongs to the FEMaster Neumann category because its solver-
 * facing contribution is entirely an RHS term. Prescribed temperatures for a
 * thermal conduction analysis are represented separately by `Temperature` in
 * the Dirichlet category.
 *
 * @see Neumann
 * @see Temperature
 * @see model::StructuralElement
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include "neumann.h"

#include <cmath>
#include <memory>
#include <string>

namespace fem::bc {

/**
 * @brief Converts nodal temperature differences to equivalent structural forces.
 *
 * `temp_field_` stores the absolute nodal temperature `T`, while `ref_temp_`
 * defines the stress-free reference temperature `T_ref`. A structural element
 * evaluates its constitutive thermal strain and corresponding equivalent nodal
 * force. For the current solid formulation, the isotropic free strain is
 *
 *     epsilon_th = alpha (T - T_ref) [1, 1, 1, 0, 0, 0]^T,
 *
 * and the element integrates the resulting thermal stress contribution through
 * a term of the form
 *
 *     f_th,e = integral_Omega_e B^T C epsilon_th dOmega.
 *
 * The exact kinematics, constitutive tangent, quadrature and sign convention
 * remain owned by the concrete structural element.
 */
struct TLoad : Neumann {
    // Types
    using Ptr = std::shared_ptr<TLoad>;

    // Scalar nodal absolute-temperature field used by the element formulations
    SPtr<model::Field> temp_field_ = nullptr;

    // Stress-free reference temperature defining Delta T = T - T_ref
    Precision ref_temp_ = NAN;

    // Construction
    TLoad() = default;
    ~TLoad() override = default;

    // Validate the nodal temperature field and delegate equivalent thermal-force
    // assembly to every structural element.
    void apply(model::ModelData& model_data, model::Field& rhs,
               Precision time, bool ignore_amplitude = false) override;

    // Diagnostics
    std::string str() const override;
};

} // namespace fem::bc
