/**
 * @file load_t.cpp
 * @brief Implements the thermal load boundary condition.
 *
 * @see src/bc/load_t.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_t.h"

#include "../core/logging.h"
#include "../model/element/element_structural.h"
#include "../model/model_data.h"

#include <sstream>

namespace fem {
namespace bc {

/**
 * @copydoc TLoad::apply
 */
void TLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time) {
    (void)time;

    logging::error(temp_field != nullptr, "Temperature field not set on TLOAD");
    logging::error(temp_field->domain == model::FieldDomain::NODE,
                   "Temperature field ", temp_field->name, " must be a node field");
    logging::error(temp_field->components == 1,
                   "Temperature field ", temp_field->name, " must have 1 component");

    for (auto& element_ptr : model_data.elements) {
        if (auto structural = element_ptr->as<model::StructuralElement>()) {
            structural->apply_tload(bc, *temp_field, ref_temp);
        }
    }
}

/**
 * @copydoc TLoad::str
 */
std::string TLoad::str() const {
    std::ostringstream os;
    os << "TLOAD: field=" << (temp_field ? temp_field->name : std::string("?"))
       << ", ref_temp=" << ref_temp;
    return os.str();
}

} // namespace bc
} // namespace fem
