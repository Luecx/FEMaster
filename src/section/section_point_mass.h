/******************************************************************************
* @file section_point_mass.h
* @brief Defines the PointMassSection class for FEM point mass sections.
*
* A point mass section includes mass, inertia, and spring properties.
*
* @see point_mass_section.cpp
* @date 21.11.2024
******************************************************************************/

#pragma once

#include "section.h"
#include "../core/types_eig.h"

namespace fem {

/******************************************************************************
* @class PointMassSection
* @brief Represents a point mass FEM section.
*
* A point mass section includes mass, rotary inertia, and spring constants.
******************************************************************************/
struct PointMassSection : Section {
   using Ptr = std::shared_ptr<PointMassSection>;

   Precision mass = 0;              ///< Mass of the point mass section.
   Vec3 rotary_inertia = Vec3::Zero(); ///< Rotary inertia for the section.
   Vec3 spring_constants = Vec3::Zero(); ///< Spring constants for translation.
   Vec3 rotary_spring_constants = Vec3::Zero(); ///< Spring constants for rotation.

   /******************************************************************************
    * @brief Outputs information about the point mass section.
    ******************************************************************************/
   void info() override;
};

} // namespace fem
