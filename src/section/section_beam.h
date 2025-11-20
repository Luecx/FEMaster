/**
* @file section_beam.h
* @brief Defines the BeamSection class for beam FEM sections.
*
* A beam section includes material, region, and additional profile information
* as well as a directional vector `n1`.
*
* @see section_beam.cpp
* @date 21.11.2024
*/

#pragma once

#include "section.h"
#include "profile.h"
#include "../core/types_eig.h"

namespace fem {

/**
* @class BeamSection
* @brief Represents a beam FEM section.
*
* A beam section includes material, region, and profile information.
*/
struct BeamSection : Section {
   using Ptr = std::shared_ptr<BeamSection>;

   Vec3 n1 = Vec3::Zero();         ///< Direction vector for the beam section (optional).
   Profile::Ptr profile = nullptr; ///< Profile associated with the beam section.

   /**
    * @brief Outputs information about the beam section.
    */
   void info() override;
};

} // namespace fem
