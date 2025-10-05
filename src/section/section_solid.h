/**
* @file section_solid.h
* @brief Defines the SolidSection class for solid FEM sections.
*
* A solid section associates a material and a region without additional
* properties.
*
* @see section_solid.cpp
* @date 21.11.2024
*/

#pragma once

#include "section.h"

namespace fem {

/**
* @class SolidSection
* @brief Represents a solid FEM section.
*
* A solid section contains material and region information.
*/
struct SolidSection : Section {
   using Ptr = std::shared_ptr<SolidSection>;

   /**
    * @brief Outputs information about the solid section.
    */
   void info() override;
};

} // namespace fem
