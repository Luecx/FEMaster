/**
 * @file load.h
 * @brief Declares the common ownership base for load-side boundary conditions.
 *
 * A `Load` stores the modifiers shared by natural and mixed boundary conditions:
 * an optional spatial coordinate system and an optional time-dependent amplitude.
 * Concrete Neumann conditions convert physical loading definitions into
 * right-hand-side contributions, while Robin conditions additionally modify the
 * system operator through their dedicated assembly interface.
 *
 * The base class deliberately retains the generic RHS-only dispatch used by
 * structural and thermal Neumann conditions. Robin conditions derive from the
 * same ownership hierarchy so heterogeneous load collectors can store them, but
 * explicitly reject the RHS-only path and require operator-assembly context.
 *
 * @see Load
 * @see Neumann
 * @see Robin
 * @see LoadCollector
 * @see amplitude.h
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "amplitude.h"
#include "bc.h"

#include "../core/printable.h"
#include "../core/types_cls.h"
#include "../core/types_num.h"
#include "../cos/coordinate_system.h"
#include "../data/field.h"

#include <memory>
#include <string>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Common state and RHS dispatch shared by load-side boundary conditions.
 *
 * Derived classes represent physical conditions that are owned by named load
 * collectors. `orientation_` defines the basis in which vector-valued load
 * components are interpreted; a null pointer means the stored components are
 * already expressed in the global model basis. `amplitude_` optionally scales
 * the nominal condition as a function of analysis time.
 *
 * Neumann conditions implement `apply()` directly and assemble only into the
 * supplied nodal right-hand-side field. Robin conditions remain in this common
 * ownership hierarchy but override the generic path with an error because their
 * temperature- or state-dependent contribution also requires access to the
 * global operator assembly.
 *
 * Every concrete load provides a compact printable representation for model
 * overviews and diagnostics.
 */
struct Load : BoundaryCondition, Printable {
    // Shared ownership type used by load collectors to store heterogeneous
    // Neumann and Robin conditions through one common base interface.
    using Ptr = std::shared_ptr<Load>;

    // Optional coordinate system in which vector-valued load components are
    // defined. A null pointer means that the stored components are already in
    // the global model basis.
    cos::CoordinateSystem::Ptr orientation_ = nullptr;

    // Optional scalar time history. A null pointer leaves the nominal load
    // unchanged; concrete implementations may also explicitly bypass scaling.
    Amplitude::Ptr amplitude_ = nullptr;

    // Enable destruction through a `Load` pointer without requiring the
    // collector to know the concrete boundary-condition type.
    virtual ~Load() = default;

    // Assemble an RHS-only boundary condition into `rhs` using model geometry
    // and topology. `time` is supplied to the optional amplitude, while
    // `ignore_amplitude` requests the unscaled nominal condition. Robin
    // conditions reject this path and require their operator-aware overload.
    virtual void apply(model::ModelData& model_data,
                       model::Field&     rhs,
                       Precision         time,
                       bool              ignore_amplitude = false) = 0;

    // Return a compact one-line description containing the concrete condition,
    // target region and parameters relevant for diagnostics.
    std::string str() const override = 0;
};

} // namespace fem::bc
