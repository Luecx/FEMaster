/**
 * @file mask_field.h
 * @brief Declares boolean selection of model-field entries.
 *
 * The matrix tools layer exposes shape-preserving masking of model fields.
 * The implementation validates the mask and returns independent field storage;
 * model::Field defines the underlying domain and value representation.
 *
 * @see mask_field
 *
 * @author Finn Eggers
 * @date 05.10.2026
 */

#pragma once

#include "../core/types_eig.h"
#include "../data/field.h"

#include <string>

namespace fem {
namespace mattools {

// Copy selected entries; unselected entries remain NaN. Preserve source shape
// and domain, require an equally shaped mask, and assign the supplied name.
model::Field mask_field(const model::Field& field,
                        const BooleanMatrix& mask,
                        const std::string& name);

} // namespace mattools
} // namespace fem
