/**
 * @file load.h
 * @brief Declares the abstract base interface for load boundary conditions.
 *
 * Concrete load types are declared in dedicated headers:
 * - `load_c.h`
 * - `load_d.h`
 * - `load_p.h`
 * - `load_v.h`
 * - `load_t.h`
 * - `load_inertial.h`
 *
 * This file only contains the polymorphic `Load` base class.
 *
 * @see src/bc/load_collector.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../core/types_cls.h"
#include "../core/types_eig.h"
#include "amplitude.h"
#include "../cos/coordinate_system.h"
#include "../data/field.h"
#include "../data/region.h"
#include "bc.h"
#include "../core/printable.h"
#include <string>

namespace fem {
namespace model {
class ModelData;
}
}

namespace fem {
namespace bc {

/**
 * @struct Load
 * @brief Base interface for load boundary conditions.
 *
 * Derived load types must implement the `apply` method to scatter their values
 * into the boundary-condition data for all entities referenced by their target
 * region.
 */
struct Load : public BoundaryCondition, public fem::Printable {
    using Ptr = std::shared_ptr<Load>; ///< Shared pointer alias for load storage.

    cos::CoordinateSystem::Ptr orientation = nullptr; ///< Optional local orientation.
    Amplitude::Ptr amplitude = nullptr;               ///< Optional time-dependent scaling.

    /**
     * @brief Virtual defaulted destructor to enable polymorphic deletion.
     */
    virtual ~Load() = default;

    /**
     * @brief Applies the load to the given model.
     *
     * @param model_data Model data that provides geometry and topology.
     * @param bc Boundary-condition storage where load contributions are added.
     */
    virtual void apply(model::ModelData& model_data, model::Field& bc, Precision time) = 0;

    /// One-line description of the load (type, target, parameters)
    std::string str() const override = 0;
};

} // namespace bc
} // namespace fem

// Convenience umbrella includes for concrete load types.
#include "load_c.h"
#include "load_d.h"
#include "load_p.h"
#include "load_v.h"
#include "load_t.h"
#include "load_inertial.h"
