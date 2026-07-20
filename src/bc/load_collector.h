/**
 * @file load_collector.h
 * @brief Declares the named collection used to construct and apply load sets.
 *
 * `LoadCollector` owns a heterogeneous sequence of `Load` objects through
 * shared pointers. It provides convenience functions that create each concrete
 * load type from parser-level parameters and a single `apply()` operation that
 * assembles the complete load set into a model field at the current analysis
 * time. The construction and dispatch logic is implemented in
 * `load_collector.cpp`.
 *
 * @see LoadCollector
 * @see Load
 * @see load_collector.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../data/collection.h"
#include "../data/region.h"
#include "load.h"

namespace fem {
namespace model {
class ModelData;
}
}

namespace fem {
namespace bc {

/**
 * @brief Stores and assembles a named heterogeneous load set.
 *
 * The inherited `model::Collection` supplies naming and storage, while this
 * class adds strongly typed factory-style insertion methods for the supported
 * boundary-condition formulations. Application remains polymorphic: each
 * stored load receives the same model, target field and analysis time and is
 * responsible for its own region traversal and integration.
 */
struct LoadCollector : model::Collection<Load::Ptr> {
    // Shared ownership type used when load sets are referenced by load cases or
    // model-level registries.
    using Ptr = std::shared_ptr<LoadCollector>;

    // Construct a named collection using the storage semantics of the generic
    // `model::Collection<Load::Ptr>` base.
    explicit LoadCollector(const std::string& name);

    // Destroy the collection and release its shared references to all stored
    // loads.
    ~LoadCollector() = default;

    // Apply every non-null load in insertion order. Contributions are
    // accumulated in `bc`, so overlapping loads naturally superimpose.
    void apply(model::ModelData& model_data, model::Field& bc, Precision time);

    // Create and store a concentrated nodal load. The generalized vector is
    // ordered as `[Fx, Fy, Fz, Mx, My, Mz]`; optional orientation and amplitude
    // objects are transferred into the new load.
    void add_cload(model::NodeRegion::Ptr region,
                   Vec6 values,
                   cos::CoordinateSystem::Ptr orientation = nullptr,
                   Amplitude::Ptr amplitude = nullptr);

    // Create and store a distributed surface traction. `values` is interpreted
    // globally unless `orientation` supplies a local basis, and `amplitude`
    // optionally scales the vector at application time.
    void add_dload(model::SurfaceRegion::Ptr region,
                   Vec3 values,
                   cos::CoordinateSystem::Ptr orientation = nullptr,
                   Amplitude::Ptr amplitude = nullptr);

    // Create and store a scalar pressure load acting normal to every surface in
    // `region`. The optional amplitude scales the nominal pressure.
    void add_pload(model::SurfaceRegion::Ptr region,
                   Precision pressure,
                   Amplitude::Ptr amplitude = nullptr);

    // Create and store a density-scaled distributed body load over structural
    // elements. An optional coordinate system may rotate the vector separately
    // at each integration point.
    void add_vload(model::ElementRegion::Ptr region,
                   Vec3 values,
                   cos::CoordinateSystem::Ptr orientation = nullptr,
                   Amplitude::Ptr amplitude = nullptr);

    // Create and store a thermal load driven by a scalar nodal temperature
    // field and a stress-free reference temperature.
    void add_tload(model::Field::Ptr temp_field, Precision ref_temp);

    // Create and store an equivalent inertia load for the rigid-body
    // acceleration field defined by translation, angular velocity and angular
    // acceleration about `center`. Point-mass features are included only when
    // explicitly requested.
    void add_inertialload(model::ElementRegion::Ptr region,
                          Vec3 center,
                          Vec3 center_acc,
                          Vec3 omega,
                          Vec3 alpha,
                          bool consider_point_masses = false);

    // Expose the stored polymorphic loads without allowing callers to replace
    // the collection itself. Individual pointed-to loads remain mutable through
    // their shared pointers.
    const std::vector<Load::Ptr>& entries() const { return this->_data; }
};
} // namespace bc
} // namespace fem
