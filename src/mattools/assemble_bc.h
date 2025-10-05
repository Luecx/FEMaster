/**
 * @file assemble_bc.h
 * @brief Provides functions for handling and assembling constraints into NodeData,
 * with options for handling duplicate entries, either summing values or enforcing
 * consistency for constraints.
 *
 * This function takes two NodeData objects and an option of how to handle duplicate
 * entries. If duplicate non-NaN values are encountered, depending on the option,
 * it either adds the entries or checks for consistent constraints. If inconsistent
 * constraints are found, it will error out.
 *
 * The function assumes that NodeData is defined as Eigen::Matrix<Precision, Eigen::Dynamic, 6>.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 * all rights reserved
 * @date Created on 28.08.2024
 */

#pragma once

#include "../core/types_eig.h"

namespace fem { namespace mattools {

/**
 * @brief Enum to define the behavior when handling duplicate entries in the
 * NodeData during the assembly process.
 */
enum class DuplicateHandling {
    ADD,    ///< Add the values of duplicate entries
    SET     ///< Set constraints and ensure consistency; throw an error if duplicate values differ
};

/**
 * @brief Assembles constraints from two NodeData objects, handling duplicate
 * entries based on the provided option.
 *
 * @param[in,out] target_bc The target NodeData into which the values from the source_bc will be added.
 * @param[in] source_bc The source NodeData that provides new values to be added or set into the target.
 * @param handling Option specifying how to handle duplicate non-NaN values (either ADD or SET).
 *
 * @throws std::runtime_error if duplicate non-NaN values are inconsistent when handling is set to SET.
 */
void assemble_bc(NodeData& target_bc, const NodeData& source_bc, DuplicateHandling handling);

} } // namespace fem::mattools
