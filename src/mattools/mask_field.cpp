/**
 * @file mask_field.cpp
 * @brief Implements boolean selection of model-field entries.
 *
 * The matrix tools layer preserves field shape and domain while copying selected
 * values into independent storage. Unselected entries remain NaN; field storage
 * and domain semantics are provided by model::Field.
 *
 * @see mask_field
 *
 * @author Finn Eggers
 * @date 05.10.2026
 */

#include "mask_field.h"

#include "../core/logging.h"

namespace fem {
namespace mattools {

/**
 * Copies selected entries into a new field with the source domain and shape.
 *
 * Shape validation precedes allocation and indexing. The destination is filled
 * with NaN so unselected entries remain undefined; selected entries retain their
 * original values, including any existing nonfinite values. The source and mask
 * are unchanged, and the result owns independent storage.
 *
 * @param field Source values and field-domain metadata.
 * @param mask Boolean selection with one entry per source row and component.
 * @param name Semantic name assigned to the result.
 * @return Independent masked field with unselected entries set to NaN.
 */
model::Field mask_field(const model::Field& field,
                        const BooleanMatrix& mask,
                        const std::string& name) {
    // Validate compatible shapes before accessing any mask entry.
    logging::error(field.rows == static_cast<Index>(mask.rows())
                   && field.components == static_cast<Index>(mask.cols()),
        "mask_field: field/mask shape mismatch");

    // Preserve metadata and initialize undefined entries before selection.
    model::Field result{
        name,
        field.domain,
        field.rows,
        field.components
    };
    result.fill_nan();

    // Copy only selected values without modifying the source field.
    for (Index row = 0; row < field.rows; ++row) {
        for (Index component = 0; component < field.components; ++component) {
            if (mask(row, component)) {
                result(row, component) = field(row, component);
            }
        }
    }

    return result;
}

} // namespace mattools
} // namespace fem
