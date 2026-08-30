/**
 * @file numerate_dofs.h
 * @brief Declares contiguous numbering of active system DOFs.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 * all rights reserved
 * @date Created on 28.08.2024
 */

#pragma once

#include "../core/types_eig.h"

namespace fem { namespace mattools {

/**
 * @brief Converts an arbitrary node-by-component activation mask into global
 * system identifiers with matching dimensions.
 *
 * Active entries are numbered contiguously from zero; inactive entries receive
 * -1. The component count is retained instead of assuming structural six-DOF
 * storage.
 *
 * @param system_dofs Boolean matrix indicating active nodal components.
 * @return Matrix of global system identifiers with the same dimensions.
 */
SystemDofIds numerate_dofs(const SystemDofs& system_dofs);

} } // namespace fem::mattools
