/**
 * @file load_collector.h
 * @brief Declares the collector for aggregating load boundary conditions.
 *
 * The `LoadCollector` extends the generic `model::Collection` to maintain a
 * polymorphic set of loads that can be applied collectively to an FEM model.
 *
 * @see src/bc/load_collector.cpp
 * @see src/bc/load.h
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
 * @struct LoadCollector
 * @brief Container that manages and applies polymorphic loads.
 *
 * By storing `Load` instances via shared pointers the collector enables users
 * to register various load types and apply all contributions in one pass.
 */
struct LoadCollector : model::Collection<Load::Ptr> {
    using Ptr = std::shared_ptr<LoadCollector>; ///< Shared pointer alias for collectors.

    /**
     * @brief Constructs a collector with the supplied name.
     *
     * @param name Identifier for the collector within the owning model.
     */
    explicit LoadCollector(const std::string& name);

    /**
     * @brief Defaulted virtual destructor for polymorphic cleanup.
     */
    ~LoadCollector() = default;

    /**
     * @brief Applies every stored load to the given model data.
     *
     * @param model_data FEM model data that provides geometry and topology.
     * @param bc Boundary-condition storage receiving all load contributions.
     */
    void apply(model::ModelData& model_data, model::Field& bc, Precision time);

    /**
     * @brief Adds a concentrated nodal load to the collector.
     *
     * @param region Node region receiving the load.
     * @param values Generalized load vector (Fx,Fy,Fz,Mx,My,Mz).
     * @param orientation Optional orientation system for interpreting the components.
     */
    void add_cload(model::NodeRegion::Ptr region,
                   Vec6 values,
                   cos::CoordinateSystem::Ptr orientation = nullptr,
                   Amplitude::Ptr amplitude = nullptr);

    /**
     * @brief Adds a distributed surface load to the collector.
     *
     * @param region Surface region receiving the traction.
     * @param values Surface traction vector.
     * @param orientation Optional orientation system for the traction components.
     */
    void add_dload(model::SurfaceRegion::Ptr region,
                   Vec3 values,
                   cos::CoordinateSystem::Ptr orientation = nullptr,
                   Amplitude::Ptr amplitude = nullptr);

    /**
     * @brief Adds a uniform pressure load to the collector.
     *
     * @param region Surface region receiving the pressure.
     * @param pressure Magnitude of the pressure.
     */
    void add_pload(model::SurfaceRegion::Ptr region,
                   Precision pressure,
                   Amplitude::Ptr amplitude = nullptr);

    /**
     * @brief Adds a volumetric load to the collector.
     *
     * @param region Element region receiving the body forces.
     * @param values Body-force components.
     * @param orientation Optional orientation system for the body-force vector.
     */
    void add_vload(model::ElementRegion::Ptr region,
                   Vec3 values,
                   cos::CoordinateSystem::Ptr orientation = nullptr,
                   Amplitude::Ptr amplitude = nullptr);

    /**
     * @brief Adds a thermal load to the collector.
     *
     * @param temp_field Temperature field that drives the thermal load.
     * @param ref_temp Reference temperature used as the unloaded state.
     */
    void add_tload(model::Field::Ptr temp_field, Precision ref_temp);

    /// Read-only access to stored loads
    const std::vector<Load::Ptr>& entries() const { return this->_data; }
};

} // namespace bc
} // namespace fem
