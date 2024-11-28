/******************************************************************************
* @file load_collector.h
* @brief Defines the LoadCollector class for managing collections of loads in an FEM model.
*
* @details The LoadCollector class extends the Collection class and provides an interface
*          to manage and apply multiple Load instances to a FEM model. Loads can be added
*          for node, element, or surface regions with specified values.
*
* @date Created on 28.11.2024
* @author Finn Eggers
******************************************************************************/

#pragma once

#include "../data/collection.h"
#include "../data/region.h"
#include "load.h"

namespace fem {

/**
* @struct LoadCollector
* @brief Manages a collection of Load objects for applying loads in an FEM model.
*
* @details The LoadCollector allows managing multiple loads for various regions (nodes,
*          elements, and surfaces) and provides functionality to apply all collected loads
*          to the FEM model. Loads are stored as shared pointers to the base Load class.
*/
struct LoadCollector : model::Collection<Load::Ptr> {
   using Ptr = std::shared_ptr<LoadCollector>;

   /**
    * @brief Constructor for the LoadCollector.
    * @param p_name Name of the load collector.
    */
   explicit LoadCollector(const std::string& p_name);

   /// Destructor
   ~LoadCollector() = default;

   /**
    * @brief Applies all loads in the collector to the FEM model.
    * @param model_data The FEM model data.
    * @param bc The boundary condition data structure.
    */
   void apply(model::ModelData& model_data, NodeData& bc);

   /**
    * @brief Adds a CLoad (nodal load) to the collector.
    * @param region Pointer to the node region.
    * @param values Load values for the degrees of freedom.
    */
   void add_cload(model::NodeRegion::Ptr region, Vec6 values);

   /**
    * @brief Adds a DLoad (distributed load) to the collector.
    * @param region Pointer to the surface region.
    * @param values Distributed load vector.
    */
   void add_dload(model::SurfaceRegion::Ptr region, Vec3 values);

   /**
    * @brief Adds a PLoad (pressure load) to the collector.
    * @param region Pointer to the surface region.
    * @param pressure Pressure value to be applied.
    */
   void add_pload(model::SurfaceRegion::Ptr region, Precision pressure);

   /**
    * @brief Adds a VLoad (volumetric load) to the collector.
    * @param region Pointer to the element region.
    * @param values Volumetric load vector.
    */
   void add_vload(model::ElementRegion::Ptr region, Vec3 values);
};

} // namespace fem
