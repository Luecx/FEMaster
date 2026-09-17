/**
 * @file load_c.h
 * @brief Defines concentrated nodal forces and moments.
 *
 * `CLoad` is the direct nodal member of the Neumann load family. It stores up to
 * three force and three moment components and adds them to every node of a
 * target region after optional amplitude scaling and coordinate transformation.
 *
 * Unlike distributed loads, no finite-element quadrature is required because
 * the physical load is already defined at discrete nodal degrees of freedom.
 *
 * @see Neumann
 * @see Load
 * @see cos::CoordinateSystem
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include "neumann.h"
#include "../../data/region.h"

#include <cmath>
#include <memory>
#include <string>

namespace fem::bc {

/**
 * @brief Applies concentrated generalized loads to a node region.
 *
 * The six nominal components follow
 *
 *     [Fx, Fy, Fz, Mx, My, Mz].
 *
 * `NaN` marks an omitted component and therefore contributes zero. With scalar
 * amplitude `a(t)` and local-to-global basis `A(x)`, a prescribed local force
 * vector is assembled as
 *
 *     f_i <- f_i + a(t) A(x_i) f_local,
 *
 * while the moment triplet is transformed by the same basis and accumulated in
 * the rotational generalized DOFs. Without an orientation, `A = I`.
 */
struct CLoad : Neumann {
    // Types
    using Ptr = std::shared_ptr<CLoad>;

    // Nominal generalized components. NaN marks an omitted force or moment
    // component and is converted to zero during assembly.
    Vec6 values_ = {NAN, NAN, NAN, NAN, NAN, NAN};

    // Target nodes receiving the same nominal generalized load
    SPtr<model::NodeRegion> region_ = nullptr;

    // Construction
    CLoad() = default;
    ~CLoad() override = default;

    // Add the transformed and amplitude-scaled generalized load to every target
    // node in the supplied RHS field.
    void apply(model::ModelData& model_data, model::Field& rhs,
               Precision time, bool ignore_amplitude = false) override;

    // Diagnostics
    std::string str() const override;
};

} // namespace fem::bc
