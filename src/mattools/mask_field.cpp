#include "mask_field.h"

#include "../core/logging.h"

namespace fem {
namespace mattools {

model::Field mask_field(const model::Field& field,
                        const BooleanMatrix& mask,
                        const std::string& name) {
    logging::error(field.rows == mask.rows()
                   && field.components == mask.cols(),
        "mask_field: field/mask shape mismatch");

    model::Field result{
        name,
        field.domain,
        field.rows,
        field.components
    };
    result.fill_nan();

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
