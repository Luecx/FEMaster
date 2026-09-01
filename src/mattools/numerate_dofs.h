/**
 * @file numerate_dofs.h
 * @brief Declares contiguous numbering of active system degrees of freedom.
 *
 * The numbering helper accepts an arbitrary node-by-component activation mask and
 * returns a system-index matrix with identical dimensions. This allows the same
 * routine to serve both the six-component structural mapping and scalar thermal
 * mappings without embedding a fixed component count in the utility interface.
 *
 * Detailed traversal and numbering semantics are documented at the definition in
 * `numerate_dofs.cpp`.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 * all rights reserved
 * @date Created on 28.08.2024
 */

#pragma once

#include "../core/types_eig.h"

namespace fem { namespace mattools {

// Convert an arbitrary node-by-component activation mask into contiguous global
// system identifiers. Active entries are numbered from zero in row-major
// traversal order; inactive entries receive -1 and the input dimensions are kept.
SystemDofIds numerate_dofs(const SystemDofs& system_dofs);

} } // namespace fem::mattools
