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

    logging::error(temp_field_ != nullptr, "Temperature field not set on TLOAD");
    logging::error(temp_field_->domain == model::FieldDomain::NODE,
                   "Temperature field ", temp_field_->name, " must be a node field");
    logging::error(temp_field_->components == 1,
                   "Temperature field ", temp_field_->name, " must have 1 component");

    // Thermal loading is element-defined: each structural element reads nodal
    // temperatures from the shared field and scatters its equivalent nodal RHS.
    for (auto& element_ptr : model_data.elements) {
        if (auto structural = element_ptr->as<model::StructuralElement>()) {
            structural->apply_tload(bc, *temp_field_, ref_temp_);
        }
    }
}

/**
 * @copydoc TLoad::str
 */
std::string TLoad::str() const {
    std::ostringstream os;

    os << "TLOAD: field="
       << (temp_field_ ? temp_field_->name : std::string("?"))
       << ", ref_temp=" << ref_temp_;

    return os.str();
}
} // namespace bc
} // namespace fem
