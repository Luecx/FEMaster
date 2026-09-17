/**
 * @file load_t.cpp
 * @brief Implements delegation of equivalent thermal-force assembly.
 *
 * A structural temperature load represents thermal expansion as an equivalent
 * mechanical RHS contribution. `TLoad` owns only the global scalar temperature
 * field and the stress-free reference temperature; each structural element owns
 * the kinematics, constitutive law and quadrature needed to convert that data to
 * nodal forces.
 *
 * For the current solid formulation, the interpolated temperature defines
 *
 *     Delta T       = T - T_ref,
 *     epsilon_th    = alpha Delta T [1, 1, 1, 0, 0, 0]^T,
 *     sigma_th      = C epsilon_th,
 *     f_th,e        = integral_Omega_e B^T sigma_th dOmega.
 *
 * Other structural element types may implement a different thermal formulation
 * or deliberately provide no thermal contribution. This dispatcher does not
 * duplicate any of that element-specific mechanics.
 *
 * @see TLoad
 * @see Neumann
 * @see model::StructuralElement
 * @see model::SolidElement
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#include "load_t.h"

#include "../../core/logging.h"
#include "../../model/element/element_structural.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>

namespace fem::bc {

/**
 * Validates the shared temperature field and delegates thermal loading to every
 * structural element.
 *
 * The source field must contain one scalar absolute temperature per node. Each
 * element gathers its own nodal values, interpolates temperature according to
 * its shape functions and constructs the corresponding equivalent nodal force.
 * The reference temperature defines the stress-free state through
 *
 *     Delta T = T - T_ref.
 *
 * The current `TLoad` representation does not multiply the temperature field by
 * the inherited scalar amplitude. Time dependence must already be represented
 * by the selected temperature field before this function is called; `time` and
 * `ignore_amplitude` are therefore intentionally unused.
 *
 * @param model_data Compiled structural elements receiving the thermal load.
 * @param rhs Generalized nodal RHS field modified by the element formulations.
 * @param time Unused by the current field-based thermal-load representation.
 * @param ignore_amplitude Unused because no additional amplitude is applied.
 */
void TLoad::apply(model::ModelData& model_data, model::Field& rhs, Precision time, bool ignore_amplitude) {
    // The prescribed nodal temperature field already represents the complete
    // thermal state for this load evaluation
    (void)time;
    (void)ignore_amplitude;

    // Validate the scalar nodal thermal field and stress-free reference state
    logging::error(temp_field_ != nullptr,
        "TLOAD: temperature field is not initialized");
    logging::error(temp_field_->domain == model::FieldDomain::NODE,
        "TLOAD: temperature field ", temp_field_->name, " must use NODE domain");
    logging::error(temp_field_->components == 1,
        "TLOAD: temperature field ", temp_field_->name, " must have one component");
    logging::error(std::isfinite(ref_temp_),
        "TLOAD: reference temperature must be finite");

    // Let each structural formulation construct its own B, constitutive thermal
    // strain and quadrature contribution. This keeps the thermal-force mapping
    // consistent with the element's actual mechanical kinematics.
    for (auto& element : model_data.elements) {
        auto* structural = element->as<model::StructuralElement>();
        if (!structural) {
            continue;
        }

        structural->apply_tload(rhs, *temp_field_, ref_temp_);
    }
}

/**
 * Builds a compact diagnostic representation of the structural thermal load.
 *
 * The output identifies the scalar nodal temperature field and the stress-free
 * reference temperature. No element interpolation or constitutive evaluation is
 * performed for diagnostics.
 *
 * @return Human-readable thermal-load description.
 */
std::string TLoad::str() const {
    std::ostringstream os;

    os << "TLOAD: field="
       << (temp_field_ ? temp_field_->name : std::string("?"))
       << ", ref_temp=" << ref_temp_;

    return os.str();
}

} // namespace fem::bc
