/**
 * @file load.h
 * @brief Declares the hierarchy of load boundary conditions.
 *
 * The BC module defines several concrete load types that can be attached to
 * node, surface, or element regions. Each load derives from `BoundaryCondition`
 * and implements the `apply` function to assemble its contribution into the
 * global boundary condition data structures.
 *
 * @see src/bc/load.cpp
 * @see src/bc/load_collector.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../core/types_cls.h"
#include "../data/node_data_dict.h"
#include "../data/region.h"
#include "bc.h"

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
struct Load : public BoundaryCondition {
    using Ptr = std::shared_ptr<Load>; ///< Shared pointer alias for load storage.

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
    virtual void apply(model::ModelData& model_data, NodeData& bc) = 0;
};

/**
 * @struct CLoad
 * @brief Concentrated nodal load.
 *
 * Applies up to six generalized force/moment components to each node in the
 * associated node region.
 */
struct CLoad : public Load {
    using Ptr = std::shared_ptr<CLoad>; ///< Shared pointer alias for concentrated loads.

    Vec6 values{NAN, NAN, NAN, NAN, NAN, NAN}; ///< Generalized load vector (Fx,Fy,Fz,Mx,My,Mz).
    model::NodeRegion::Ptr region = nullptr;   ///< Target node region.

    /**
     * @brief Default constructor.
     */
    CLoad() = default;

    /**
     * @brief Defaulted virtual destructor for polymorphic cleanup.
     */
    ~CLoad() override = default;

    /**
     * @brief Applies the load to all nodes inside `region`.
     *
     * @param model_data Provides access to nodal topology (unused).
     * @param bc Node-based boundary-condition storage.
     */
    void apply(model::ModelData& model_data, NodeData& bc) override;
};

/**
 * @struct DLoad
 * @brief Distributed surface load with vector traction.
 *
 * Applies translational forces to the nodes of each surface patch referenced by
 * the associated region.
 */
struct DLoad : public Load {
    using Ptr = std::shared_ptr<DLoad>; ///< Shared pointer alias for distributed loads.

    Vec3 values{NAN, NAN, NAN};                ///< Surface traction components.
    model::SurfaceRegion::Ptr region = nullptr; ///< Target surface region.

    /**
     * @brief Default constructor.
     */
    DLoad() = default;

    /**
     * @brief Defaulted virtual destructor for polymorphic cleanup.
     */
    ~DLoad() override = default;

    /**
     * @brief Applies the distributed traction to all surfaces inside `region`.
     *
     * @param model_data Provides surface connectivity and geometry.
     * @param bc Node-based boundary-condition storage.
     */
    void apply(model::ModelData& model_data, NodeData& bc) override;
};

/**
 * @struct PLoad
 * @brief Uniform surface pressure load.
 *
 * Distributes a scalar pressure to the nodes of each surface in the associated
 * region.
 */
struct PLoad : public Load {
    using Ptr = std::shared_ptr<PLoad>; ///< Shared pointer alias for pressure loads.

    Precision pressure{NAN};                  ///< Magnitude of the applied pressure.
    model::SurfaceRegion::Ptr region = nullptr; ///< Target surface region.

    /**
     * @brief Default constructor.
     */
    PLoad() = default;

    /**
     * @brief Defaulted virtual destructor for polymorphic cleanup.
     */
    ~PLoad() override = default;

    /**
     * @brief Applies the pressure to all surfaces inside `region`.
     *
     * @param model_data Provides surface connectivity and geometry.
     * @param bc Node-based boundary-condition storage.
     */
    void apply(model::ModelData& model_data, NodeData& bc) override;
};

/**
 * @struct VLoad
 * @brief Volumetric load applied to structural elements.
 *
 * Applies a body-force vector to each structural element in the associated
 * element region.
 */
struct VLoad : public Load {
    using Ptr = std::shared_ptr<VLoad>; ///< Shared pointer alias for volumetric loads.

    Vec3 values{NAN, NAN, NAN};                 ///< Body-force components.
    model::ElementRegion::Ptr region = nullptr; ///< Target element region.

    /**
     * @brief Default constructor.
     */
    VLoad() = default;

    /**
     * @brief Defaulted virtual destructor for polymorphic cleanup.
     */
    ~VLoad() override = default;

    /**
     * @brief Applies the body forces to all elements inside `region`.
     *
     * @param model_data Provides access to the element container.
     * @param bc Node-based boundary-condition storage.
     */
    void apply(model::ModelData& model_data, NodeData& bc) override;
};

/**
 * @struct TLoad
 * @brief Thermal load derived from a temperature field.
 *
 * Generates equivalent nodal forces by feeding a temperature field to each
 * structural element.
 */
struct TLoad : public Load {
    using Ptr = std::shared_ptr<TLoad>; ///< Shared pointer alias for thermal loads.

    model::NodeField::Ptr temp_field = nullptr; ///< Temperature field reference.
    Precision ref_temp{NAN};                    ///< Reference temperature for zero load.

    /**
     * @brief Default constructor.
     */
    TLoad() = default;

    /**
     * @brief Defaulted virtual destructor for polymorphic cleanup.
     */
    ~TLoad() override = default;

    /**
     * @brief Applies the thermal load to all structural elements in the model.
     *
     * @param model_data Provides access to the element container.
     * @param bc Node-based boundary-condition storage.
     */
    void apply(model::ModelData& model_data, NodeData& bc) override;
};

} // namespace bc
} // namespace fem
