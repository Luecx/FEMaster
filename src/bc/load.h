/**
 * @file load.h
 * @brief Declares the polymorphic interface shared by structural load types.
 *
 * Structural loads are Neumann boundary data that assemble generalized nodal
 * force contributions. Thermal Neumann conditions use their own scalar assembly
 * interface under `bc/neumann` because convection contributes to both matrix and
 * right-hand side.
 *
 * @see Neumann
 * @see LoadCollector
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "amplitude.h"
#include "neumann/neumann.h"

#include "../core/types_cls.h"
#include "../core/types_eig.h"
#include "../cos/coordinate_system.h"
#include "../data/field.h"
#include "../data/region.h"

#include <string>

namespace fem {
namespace model {
struct ModelData;
}
}

namespace fem {
namespace bc {

/**
 * @brief Defines the common assembly interface for structural Neumann loads.
 */
struct Load : public Neumann {
    using Ptr = std::shared_ptr<Load>;

    cos::CoordinateSystem::Ptr orientation_ = nullptr;
    Amplitude::Ptr             amplitude_   = nullptr;

    ~Load() override = default;

    virtual void apply(model::ModelData& model_data,
                       model::Field&     bc,
                       Precision         time,
                       bool              ignore_amplitude = false) = 0;

    std::string str() const override = 0;
};
} // namespace bc
} // namespace fem

#include "load_c.h"
#include "load_d.h"
#include "load_p.h"
#include "load_v.h"
#include "load_t.h"
#include "load_inertial.h"
