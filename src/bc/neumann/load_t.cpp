/**
 * @file load_t.cpp
 * @brief Implements structural thermal-expansion loading.
 *
 * Thermal-force generation depends on each structural element's kinematics,
 * material law and interpolation. This implementation therefore validates the
 * shared nodal temperature field and forwards it, together with the stress-free
 * reference temperature, to every structural element for element-specific
 * equivalent-force assembly.
 *
 * @see load_t.h
 * @see model::StructuralElement
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_t.h"

#include "../../core/logging.h"
#include "../../model/element/element_structural.h"
#include "../../model/model_data.h"

#include <sstream>

namespace fem::bc {

/**
 * Forwards the prescribed nodal temperature field to structural elements.
 *
 * Each structural element computes its own equivalent thermal force from
 * interpolation, kinematics, material behavior and the stress-free reference
 * temperature. The common time and amplitude arguments are unused because the
 * current temperature state is represented completely by `temp_field_`.
 *
 * @param model_data Model topology and structural elements.
 * @param rhs Generalized nodal right-hand-side field receiving thermal force.
 * @param time Unused analysis time retained by the common load interface.
 * @param ignore_amplitude Unused common-interface amplitude flag.
 */
void TLoad::apply(model::ModelData& model_data,
                  model::Field&     rhs,
                  Precision         time,
                  bool              ignore_amplitude) {
    // The current structural thermal-load representation receives its complete
    // state from `temp_field_`; no additional scalar time or amplitude scaling is
    // applied in this implementation.
    (void) time;
    (void) ignore_amplitude;

    // Validate the scalar nodal temperature field and structural target layout
    // before any element attempts to gather temperatures or assemble forces.
    logging::error(temp_field_ != nullptr,
        "Temperature field not set on TLOAD");
    logging::error(temp_field_->domain == model::FieldDomain::NODE,
        "Temperature field ", temp_field_->name, " must be a node field");
    logging::error(temp_field_->components == 1,
        "Temperature field ", temp_field_->name, " must have 1 component");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components >= 6,
        "TLoad: target field must be a NODE field with at least six components");

    // Traverse the complete compiled element container because TLOAD is defined
    // by its temperature field rather than by a separate element region. Only
    // structural elements consume thermal expansion and generate mechanical RHS
    // contributions; all other element capabilities are ignored.
    for (auto& element : model_data.elements) {
        if (!element) continue;
        if (auto structural = element->as<model::StructuralElement>()) {
            structural->apply_tload(rhs, *temp_field_, ref_temp_);
        }
    }
}

/**
 * Builds the diagnostic representation of the structural thermal load.
 *
 * The result identifies the source temperature field and reports the stress-free
 * reference temperature used by structural element evaluation.
 *
 * @return Human-readable load description.
 */
std::string TLoad::str() const {
    std::ostringstream os;

    // Report the source field by name when available and always include the
    // reference temperature defining the stress-free configuration.
    os << "TLOAD: field="
       << (temp_field_ ? temp_field_->name : std::string("?"))
       << ", ref_temp=" << ref_temp_;

    return os.str();
}

} // namespace fem::bc
