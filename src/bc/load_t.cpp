/**
 * @file load_t.cpp
 * @brief Implements delegation of nodal temperature loading to structural elements.
 *
 * Thermal-force generation depends on each element's kinematics, material law
 * and interpolation. This implementation therefore validates the shared nodal
 * temperature field and forwards it, together with the stress-free reference
 * temperature, to every structural element for element-specific assembly.
 *
 * @see load_t.h
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
 * Forwards the prescribed nodal temperature field to structural elements.
 *
 * Each element computes its own equivalent thermal force from interpolation,
 * kinematics, material behavior and the stress-free reference temperature.
 * The common time and amplitude arguments are unused by this representation.
 *
 * @param model_data Model topology and structural elements.
 * @param bc Generalized nodal field receiving thermal force.
 * @param time Unused analysis time retained by the common interface.
 * @param ignore_amplitude Unused common-interface flag.
 */
void TLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude) {
    // The current thermal-load representation receives its complete state from
    // `temp_field_`; neither an additional scalar time nor amplitude scaling is
    // applied in this implementation.
    (void)time;
    (void)ignore_amplitude;

    // A thermal load requires a scalar temperature value at every relevant node.
    // Validate the pointer, field domain and component count before any element
    // attempts to gather values from it.
    logging::error(temp_field_ != nullptr, "Temperature field not set on TLOAD");
    logging::error(temp_field_->domain == model::FieldDomain::NODE,
                   "Temperature field ", temp_field_->name, " must be a node field");
    logging::error(temp_field_->components == 1,
                   "Temperature field ", temp_field_->name, " must have 1 component");

    // Traverse the complete element container because thermal loading is not
    // restricted by a separate region in `TLoad`. Non-structural elements are
    // ignored; structural elements gather their own nodal temperatures and use
    // their section and material data to form the thermal contribution.
    for (auto& element_ptr : model_data.elements) {
        if (auto structural = element_ptr->as<model::StructuralElement>()) {
            structural->apply_tload(bc, *temp_field_, ref_temp_);
        }
    }
}

/**
 * Builds the diagnostic representation of the thermal load.
 *
 * The result identifies the source temperature field and reports the stress-free
 * reference temperature used by structural element evaluation.
 *
 * @return Human-readable load description.
 */
std::string TLoad::str() const {
    std::ostringstream os;

    // Report the source field by name when available and always include the
    // reference temperature that defines the stress-free configuration.
    os << "TLOAD: field="
       << (temp_field_ ? temp_field_->name : std::string("?"))
       << ", ref_temp=" << ref_temp_;

    return os.str();
}
} // namespace bc
} // namespace fem
