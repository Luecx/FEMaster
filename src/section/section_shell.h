/**
* @file section_shell.h
* @brief Defines the ShellSection class for shell FEM sections.
*
* A shell section includes material, region, and thickness information.
*
* @see shell_section.cpp
* @date 21.11.2024
*/

#pragma once

#include "section.h"
#include "../core/types_eig.h"

namespace fem {

/**
* @class ShellSection
* @brief Represents a shell FEM section.
*
* A shell section includes material, region, and thickness properties.
*/
struct ShellSection : Section {
   using Ptr = std::shared_ptr<ShellSection>;

   Precision thickness = 1.0; ///< Thickness of the shell section.

   /**
    * @brief Outputs information about the shell section.
    */
   void info() override;

   /**
    * @brief One-line section description.
    */
    std::string str() const override;
};

} // namespace fem
