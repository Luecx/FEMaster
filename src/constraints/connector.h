/**
 * @file connector.h
 * @brief Declares connector constraints that tie degrees of freedom.
 *
 * Connectors relate the motion of two nodes through a user-defined local
 * coordinate system and a bit-mask describing which generalized DOFs are
 * constrained.
 *
 * @see src/constraints/connector.cpp
 * @see src/constraints/equation.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../core/types_cls.h"
#include "../core/types_eig.h"
#include "../cos/coordinate_system.h"
#include "equation.h"

namespace fem {
namespace model {
class ModelData;
}
}

namespace fem {
namespace constraint {

/**
 * @enum ConnectorType
 * @brief Bit mask that defines which generalized DOFs are constrained.
 *
 * Bit layout (LSB to MSB): Tx, Ty, Tz, Rx, Ry, Rz.
 */
enum ConnectorType : int {
    None        = 0b000000, ///< No degree of freedom is constrained.
    Beam        = 0b111111, ///< All six DOFs are constrained.
    Hinge       = 0b111011, ///< Rotation about the local x-axis remains free.
    Cylindrical = 0b011011, ///< Translation along and rotation about local x-axis are free.
    Translator  = 0b011111, ///< Only translation in local x-direction remains free.
    Join        = 0b111000, ///< Rotations remain free, translations constrained.
    JoinRx      = 0b111100  ///< Rotation about local x-axis remains free in addition to translations.
};

/**
 * @class Connector
 * @brief Couples two nodes according to a connector type and local frame.
 */
class Connector {
    fem::cos::CoordinateSystem::Ptr coordinate_system_; ///< Local coordinate system for projection.
    ConnectorType connector_type_; ///< Connector type describing constrained DOFs.
    Dofs constrained_dofs_;        ///< Per-DOF constraint mask.
    ID node_1_id_;                 ///< First node identifier.
    ID node_2_id_;                 ///< Second node identifier.

public:
    /**
     * @brief Constructs a connector for a node pair.
     *
     * @param node_1_id Identifier of the first node.
     * @param node_2_id Identifier of the second node.
     * @param coordinate_system Local frame used to interpret the connector type.
     * @param connector_type Bit mask describing constrained DOFs.
     */
    Connector(ID node_1_id,
              ID node_2_id,
              fem::cos::CoordinateSystem::Ptr coordinate_system,
              ConnectorType connector_type);

    /**
     * @brief Virtual defaulted destructor.
     */
    virtual ~Connector() = default;

    /**
     * @brief Generates constraint equations produced by this connector.
     *
     * @param system_nodal_dofs Global DOF indices per node.
     * @param model_data Reference to model data (coordinates are required).
     * @return Equations Collection of linear constraint equations.
     */
    Equations get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data);

    /**
     * @brief Computes the impacted global DOFs considering the local frame.
     *
     * @return Dofs Mask describing affected generalized DOFs.
     */
    [[nodiscard]] Dofs dofs() const;

    /**
     * @brief Provides the first node identifier.
     *
     * @return ID ID of the first node.
     */
    [[nodiscard]] ID node_1() const { return node_1_id_; }

    /**
     * @brief Provides the second node identifier.
     *
     * @return ID ID of the second node.
     */
    [[nodiscard]] ID node_2() const { return node_2_id_; }
};

} // namespace constraint
} // namespace fem
