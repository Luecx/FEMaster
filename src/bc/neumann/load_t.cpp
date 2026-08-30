/**
 * @file load_t.cpp
 * @brief Implements delegation of nodal temperature loading to structural elements.
 *
 * Thermal-force generation depends on each element's kinematics, material law
 * and interpolation. This implementation validates the shared nodal temperature
 * field and forwards it with the stress-free reference temperature.
 *
 * The generated equivalent force modifies only the structural right-hand side.
 * The optional system DOF map and LHS triplet list of the common load interface
 * are therefore intentionally unused.
 *
 * @see load_t.h
 * @author Finn Eggers
 * @date 30.08.2026
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
 * @param model_data Model topology and structural elements.
 * @param rhs Structural nodal RHS receiving equivalent thermal force.
 * @param time Unused analysis time retained by the common load interface.
 * @param ignore_amplitude Unused common-interface amplitude flag.
 * @param system_dof_ids Unused optional system DOF map.
 * @param lhs Unused optional system-matrix triplet list.
 */
void TLoad::apply(model::ModelData&       model_data,
                  model::Field&           rhs,
                  Precision               time,
                  bool                    ignore_amplitude,
                  const SystemDofIds*      system_dof_ids,
                  TripletList*             lhs) {
    (void) time;
    (void) ignore_amplitude;
    (void) system_dof_ids;
    (void) lhs;

    logging::error(temp_field_ != nullptr,
        "Temperature field not set on TLOAD");
    logging::error(temp_field_->domain == model::FieldDomain::NODE,
        "Temperature field ", temp_field_->name, " must be a node field");
    logging::error(temp_field_->components == 1,
        "Temperature field ", temp_field_->name, " must have 1 component");

    // Let each structural element assemble its formulation-specific thermal force.
    for (auto& element : model_data.elements) {
        if (auto structural = element->as<model::StructuralElement>()) {
            structural->apply_tload(rhs, *temp_field_, ref_temp_);
        }
    }
}

/**
 * Builds the diagnostic representation of the structural thermal load.
 *
 * @return Human-readable load description.
 */
std::string TLoad::str() const {
    std::ostringstream os;
    os << "TLOAD: field="
       << (temp_field_ ? temp_field_->name : std::string("?"))
       << ", ref_temp=" << ref_temp_;
    return os.str();
}

} // namespace fem::bc
