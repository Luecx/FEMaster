/******************************************************************************
* @file version.h
* @brief Defines the version number for the FEMaster solver.
*
* This header file contains the version information (major, minor, patch)
* for the current release of the FEMaster solver. This versioning follows
* Semantic Versioning principles, where:
* - MAJOR version increments signify incompatible API changes,
* - MINOR version increments add backward-compatible functionality, and
* - PATCH version increments apply backward-compatible bug fixes.
*
* @author Finn Eggers
* @date 24.10.2024
******************************************************************************/

#pragma once

namespace fem {

/******************************************************************************
* @brief The major version number for the FEMaster solver.
*
* This number increments when incompatible API changes are introduced.
******************************************************************************/
#define VERSION_MAJOR 1

/******************************************************************************
* @brief The minor version number for the FEMaster solver.
*
* This number increments when backward-compatible functionality is added.
******************************************************************************/
#define VERSION_MINOR 0

/******************************************************************************
* @brief The patch version number for the FEMaster solver.
*
* This number increments for backward-compatible bug fixes.
******************************************************************************/
#define VERSION_PATCH 1

} // namespace fem

