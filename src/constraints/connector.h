/******************************************************************************
* @file connector.h
* @brief Defines the Connector class for applying constraints between nodes using a connector element.
*
* The Connector class links two nodes and constrains specific degrees of freedom (DOFs) using a
* custom coordinate system. The constrained DOFs are defined by the connector type, which can
* represent various physical connections such as beams or hinges.
*
* @details This class supports _connectors that constrain translations and/or rotations in specific
* directions based on the connector type, which is defined using a bitmask.
*
* @see connector.cpp
* @author Finn Eggers
* @date 08.10.2024
******************************************************************************/

#pragma once  // Ensures this file is only included once during compilation

#include "../cos/coordinate_system.h"  // For custom coordinate system handling
#include "../core/types_eig.h"             // For general types (ID, Triplet, etc.)

#include "equation.h"

namespace fem::constraint {

/******************************************************************************
* @enum ConnectorType
* @brief Enum for defining different types of _connectors using bitwise integer representations.
*
* Each connector type is represented using a bitmask, where each bit represents a degree of freedom (DOF):
* - Bit 0: Translation-X
* - Bit 1: Translation-Y
* - Bit 2: Translation-Z
* - Bit 3: Rotation-X
* - Bit 4: Rotation-Y
* - Bit 5: Rotation-Z
******************************************************************************/
enum ConnectorType : int {
   None        = 0b000000,    ///< No DOFs are constrained.
   Beam        = 0b111111,    ///< All 6 DOFs are constrained.
   Hinge       = 0b111011,    ///< All DOFs except rotation around X (bit 3 is 0).
   Cylindrical = 0b011011,    ///< Translation and rotation around X are allowed.
   Translator  = 0b011111,    ///< Only translation in the X direction is allowed.
};

/******************************************************************************
* @class Connector
* @brief Class for representing a connector element that links two nodes.
*
* The Connector class constrains specific degrees of freedom (DOFs) between two nodes
* using a custom coordinate system. The DOFs to be constrained are defined by the
* connector type, which determines the physical connection behavior (e.g., beam or hinge).
******************************************************************************/
class Connector {
   fem::cos::CoordinateSystem::Ptr coordinate_system_;    ///< Pointer to the custom coordinate system for transforming DOFs.
   ConnectorType                   connector_type_;       ///< Type of the connector (e.g., Beam, Hinge).
   Dofs                            constrained_dofs_;     ///< Vector indicating which DOFs are constrained (6 DOFs).
   ID                              node_1_id_;            ///< ID of the first node.
   ID                              node_2_id_;            ///< ID of the second node.

   public:
   /******************************************************************************
    * @brief Constructor for the Connector class.
    *
    * Initializes a Connector object that links two nodes and constrains the specified
    * degrees of freedom (DOFs) based on the connector type. The coordinate system
    * is used to transform the local DOF directions to the global system.
    *
    * @param node_1_id ID of the first node.
    * @param node_2_id ID of the second node.
    * @param coordinate_system Pointer to the custom coordinate system for the connector.
    * @param connector_type The type of the connector, which defines the constrained DOFs.
    ******************************************************************************/
   Connector(ID node_1_id,
             ID node_2_id,
             fem::cos::CoordinateSystem::Ptr coordinate_system,
             ConnectorType connector_type);

   virtual ~Connector() = default;    ///< Default destructor for the Connector class.

   /******************************************************************************
    * @brief Generates the constraint equations for the connector element.
    *
    * This method generates the constraint equations between the two nodes based
    * on the specified constrained degrees of freedom (DOFs) and their global
    * system DOF indices. The equations are transformed using the custom coordinate
    * system for accurate coupling in the global system.
    *
    * @param system_nodal_dofs The global system DOF IDs for each node.
    * @param node_coords The coordinates of all nodes (unused in this implementation).
    * @param row_offset The starting row for inserting the constraint equations.
    * @return TripletList A list of triplets representing the constraint equations.
    ******************************************************************************/
    Equations get_equations(SystemDofIds& system_nodal_dofs, model::ModelData& model_data);

};

}    // namespace fem::constraint
