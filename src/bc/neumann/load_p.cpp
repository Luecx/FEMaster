/**
 * @file load_p.cpp
 * @brief Implements pressure integration along geometric surface normals.
 *
 * Each selected surface converts the scalar pressure into a vector traction at
 * its quadrature points, using its geometric normal, and distributes the
 * integrated contribution consistently to its nodes.
 *
 * Pressure contributes only to the structural right-hand side. The optional DOF
 * map and LHS triplet list of the common load interface are intentionally unused.
 *
 * @see load_p.h
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "load_p.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <sstream>

namespace fem::bc {

/**
 * Integrates scalar pressure over the selected surfaces.
 *
 * @param model_data Global nodal positions required for surface geometry.
 * @param rhs Structural nodal RHS receiving the consistent pressure force.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether amplitude scaling is disabled.
 * @param system_dof_ids Unused optional system DOF map.
 * @param lhs Unused optional system-matrix triplet list.
 */
void PLoad::apply(model::ModelData&       model_data,
                  model::Field&           rhs,
                  Precision               time,
                  bool                    ignore_amplitude,
                  const SystemDofIds*      system_dof_ids,
                  TripletList*             lhs) {
    (void) system_dof_ids;
    (void) lhs;

    logging::error(model_data.positions != nullptr,
        "positions field not set in model data");
    logging::error(region_ != nullptr,
        "PLoad: target surface region not set");
    const auto& node_positions = *model_data.positions;

    const Precision scale           = amplitude_ && !ignore_amplitude ? amplitude_->evaluate(time) : Precision(1);
    const Precision scaled_pressure = pressure_ * scale;

    for (const ID surf_id : *region_) {
        auto surface = model_data.surfaces[surf_id];
        if (!surface) {
            continue;
        }

        // Surface quadrature supplies the physical normal and consistent nodal weighting.
        surface->integrate_vector_field(
            node_positions,
            rhs,
            [&](const Vec3& position) -> Vec3 {
                const Vec2 local = surface->global_to_local(position, node_positions);
                return -scaled_pressure * surface->normal(node_positions, local);
            }
        );
    }
}

/**
 * Builds the diagnostic representation of the pressure load.
 *
 * @return Human-readable load description.
 */
std::string PLoad::str() const {
    std::ostringstream os;
    os << "PLOAD: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", p=" << pressure_;

    if (amplitude_)
        os << ", amplitude=" << amplitude_->name;

    return os.str();
}

} // namespace fem::bc
