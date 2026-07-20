/**
 * @file load.h
 * @brief Declares the polymorphic interface shared by all load types.
 *
 * A `Load` converts a physical loading definition into contributions to a
 * model field. Concrete implementations target nodes, surfaces, elements or
 * all structural elements and may optionally interpret their components in a
 * spatial coordinate system or scale them through a time-dependent amplitude.
 * The individual load formulations are declared in their dedicated headers,
 * which are included at the end of this umbrella header for convenient access.
 *
 * @see Load
 * @see LoadCollector
 * @see amplitude.h
 * @see load_collector.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "amplitude.h"
#include "bc.h"

#include "../core/printable.h"
#include "../core/types_cls.h"
#include "../core/types_eig.h"
#include "../cos/coordinate_system.h"
#include "../data/field.h"
#include "../data/region.h"

#include <string>

namespace fem {
namespace model {
class ModelData;
}
}

namespace fem {
namespace bc {

/**
 * @brief Defines the common assembly interface for physical load conditions.
 *
 * Derived classes implement `apply()` to integrate or scatter their physical
 * load into the supplied boundary-condition field. The base class stores the
 * two optional modifiers common to several formulations: a coordinate system
 * for interpreting vector components and an amplitude for temporal scaling.
 * Each load also provides a compact printable representation for diagnostics
 * and model summaries.
 */
struct Load : public BoundaryCondition, public fem::Printable {
    // Shared ownership type used by `LoadCollector` to hold heterogeneous load
    // implementations in one contiguous collection of pointers.
    using Ptr = std::shared_ptr<Load>;

    // Optional coordinate system in which vector-valued load components are
    // defined. A null pointer means that the stored components are already in
    // the global model basis.
    cos::CoordinateSystem::Ptr orientation_ = nullptr;

    // Optional scalar time history. A null pointer leaves the nominal load
    // unchanged; concrete implementations may also explicitly bypass scaling.
    Amplitude::Ptr amplitude_ = nullptr;

    // Enable destruction through a `Load` pointer without requiring each
    // collector to know the concrete load type.
    virtual ~Load() = default;

    // Assemble this load into `bc` using geometry and topology from
    // `model_data`. `time` is passed to the optional amplitude, while
    // `ignore_amplitude` requests assembly of the unscaled nominal load.
    virtual void apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude = false) = 0;

    // Return a compact one-line description containing the concrete load type,
    // target region and the parameters relevant for diagnostics.
    std::string str() const override = 0;
};
} // namespace bc
} // namespace fem

// Expose all concrete load declarations through the common load header. The
// concrete headers include this file themselves, relying on `#pragma once` to
// terminate the cyclic umbrella inclusion after the base interface is known.
#include "load_c.h"
#include "load_d.h"
#include "load_p.h"
#include "load_v.h"
#include "load_t.h"
#include "load_inertial.h"
