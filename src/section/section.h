/**
* @file section.h
* @brief Base class definition for different FEM section types.
*
* This file contains the abstract base class Section and its derived classes,
* including SolidSection, BeamSection, ShellSection, and PointMassSection.
* Each derived class represents a specific type of section with associated
* properties and behavior.
*
* @see section.cpp
* @see solid_section.h
* @see beam_section.h
* @see shell_section.h
* @see point_mass_section.h
*
* @author Finn Eggers
* @date 21.11.2024
*/

#pragma once

#include "../material/material.h"
#include "../data/region.h"
#include <memory>

namespace fem {

/**
* @class Section
* @brief Abstract base class for FEM sections.
*
* The Section class provides a common interface and base implementation
* for different section types in the FEM solver. Each derived class
* specializes the section properties and behavior.
*/
struct Section {
   using Ptr = std::shared_ptr<Section>;

   virtual ~Section() = default; ///< Default virtual destructor for polymorphism.

   material::Material::Ptr material = nullptr; ///< Material associated with the section.
   model::ElementRegion::Ptr region = nullptr; ///< Element region associated with the section.

   /**
    * @brief Dynamic cast helper method.
    *
    * Casts this Section to a derived type.
    *
    * @tparam T Derived type to cast to.
    * @return T* Pointer to the derived type or nullptr if the cast fails.
    */
   template<typename T>
   T* as() {
       return dynamic_cast<T*>(this);
   }

   /**
    * @brief Outputs information about the section.
    */
   virtual void info() = 0;
};

} // namespace fem
