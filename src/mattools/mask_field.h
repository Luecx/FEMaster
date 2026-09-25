#pragma once

#include "../core/types_eig.h"
#include "../data/field.h"

#include <string>

namespace fem {
namespace mattools {

/**
 * Returns a copy-shaped field containing only entries selected by a boolean mask.
 *
 * Entries for which mask(row, component) is false are set to NaN. The mask must
 * have the same row/component shape as the source field. The field domain is
 * preserved while the caller supplies the semantic name of the result.
 */
model::Field mask_field(const model::Field& field,
                        const BooleanMatrix& mask,
                        const std::string& name);

} // namespace mattools
} // namespace fem
