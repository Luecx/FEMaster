/**
 * @file load_inertial.cpp
 * @brief Implements equivalent forces for prescribed rigid-body acceleration.
 *
 * Distributed structural mass is handled through each element's density-scaled
 * vector-field integrator. Optional `PointMass` features are processed
 * separately because their translational mass and diagonal rotary inertia are
 * concentrated at nodes rather than represented by structural quadrature.
 *
 * @see load_inertial.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_inertial.h"

#include "../core/logging.h"
#include "../feature/point_mass.h"
#include "../model/element/element_structural.h"
#include "../model/model_data.h"

#include <Eigen/Geometry>

#include <sstream>

namespace fem {
namespace bc {
namespace {

/**
 * Adds equivalent rigid-body inertia forces and moments for point masses.
 *
 * Translational acceleration includes reference-point, tangential and
 * centripetal terms. Translational mass and diagonal rotary inertia are then
 * assembled directly into the associated node's six generalized components.
 *
 * @param model_data Point-mass features and nodal position data.
 * @param bc Generalized nodal field receiving the result.
 * @param center Reference point of the rigid-body motion.
 * @param center_acc Translational acceleration at the reference point.
 * @param omega Angular velocity.
 * @param alpha Angular acceleration.
 */
void apply_point_mass_inertial_contribution(model::ModelData& model_data,
                                            model::Field& bc,
                                            const Vec3& center,
                                            const Vec3& center_acc,
                                            const Vec3& omega,
                                            const Vec3& alpha) {
    // Point-mass translational acceleration depends on the global nodal position
    // relative to the prescribed rigid-body reference center.
    logging::error(model_data.positions != nullptr,
                   "InertialLoad (point masses): positions field not set in model data");
    const auto& positions = *model_data.positions;

    // Search all generic model features and retain only actual point masses with
    // a valid target node region.
    for (const auto& feature_ptr : model_data.features) {
        const auto* pm = dynamic_cast<const feature::PointMass*>(feature_ptr.get());
        if (!pm || !pm->region_) {
            continue;
        }

        // A point-mass feature may assign the same mass and rotary inertia to
        // multiple nodes through its region.
        for (const ID node_id : *pm->region_) {
            // Validate the external identifier before converting it to the index
            // type expected by the position and load fields.
            logging::error(node_id >= 0 && static_cast<Index>(node_id) < positions.rows,
                           "InertialLoad (point masses): node id ", node_id,
                           " out of bounds for positions field with ", positions.rows, " rows");

            const Index node = static_cast<Index>(node_id);
            const Vec3 x = positions.row_vec3(node);
            const Vec3 r = x - center;

            // Rigid-body acceleration at the node consists of the translational
            // center acceleration, the tangential term `alpha x r` and the
            // centripetal term `omega x (omega x r)`. No Coriolis term appears
            // because the mass point has no prescribed relative velocity here.
            const Vec3 a_rot = alpha.cross(r) + omega.cross(omega.cross(r));

            // D'Alembert's equivalent inertia force opposes the prescribed
            // acceleration. Accumulate rather than assign so multiple loads and
            // features superimpose in the same nodal field.
            const Vec3 dF = -pm->mass_ * (center_acc + a_rot);
            bc(node, 0) += dF(0);
            bc(node, 1) += dF(1);
            bc(node, 2) += dF(2);

            // Rotational terms can be stored only when the target field exposes
            // the three moment columns in addition to translational forces.
            if (bc.components >= 6) {
                // `rotary_inertia_` stores diagonal principal inertia values in
                // the basis used by the prescribed angular vectors. Component-
                // wise products evaluate J*alpha and J*omega without forming a
                // dense inertia tensor.
                const Vec3 Jalpha = pm->rotary_inertia_.cwiseProduct(alpha);
                const Vec3 Jomega = pm->rotary_inertia_.cwiseProduct(omega);

                // Euler's rigid-body moment is J*alpha + omega x (J*omega). Its
                // negative is assembled as the equivalent inertia moment.
                const Vec3 dM = -(Jalpha + omega.cross(Jomega));
                bc(node, 3) += dM(0);
                bc(node, 4) += dM(1);
                bc(node, 5) += dM(2);
            }
        }
    }
}
} // namespace

/**
 * Assembles equivalent inertia loads for structural mass and point masses.
 *
 * Structural elements integrate density-scaled rigid-body acceleration over
 * their volume. Point-mass features are handled separately because their mass
 * and rotary inertia are concentrated at nodes.
 *
 * @param model_data Model topology, fields and mass-related features.
 * @param bc Generalized nodal field receiving the contribution.
 * @param time Unused analysis time retained by the common interface.
 * @param ignore_amplitude Unused common-interface flag.
 */
void InertialLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude) {
    // The current inertial definition is time-independent at this level and does
    // not use the amplitude inherited from `Load`.
    (void)time;
    (void)ignore_amplitude;
    logging::error(region_ != nullptr, "InertialLoad: region not set");

    // Structural elements in the selected region integrate their distributed
    // mass against the prescribed acceleration field.
    for (const ID el_id : *region_) {
        auto& el_ptr = model_data.elements[el_id];
        if (!el_ptr) {
            continue;
        }

        // Non-structural entries do not provide mass-weighted vector-field
        // integration and therefore cannot contribute here.
        auto structural = el_ptr->as<model::StructuralElement>();
        if (!structural) {
            continue;
        }

        // Capture the rigid-body parameters by value so the callback is
        // self-contained while the element evaluates it at its integration
        // points. The returned vector is acceleration with the D'Alembert sign;
        // density scaling is supplied by the element integrator.
        auto f = [c = center_, a0 = center_acc_, w = omega_, al = alpha_](const Vec3& x) -> Vec3 {
            const Vec3 r = x - c;
            const Vec3 a_rot = al.cross(r) + w.cross(w.cross(r));
            return -(a0 + a_rot);
        };
        structural->integrate_vector_field(bc, /*scale_by_density=*/true, f);
    }

    // Point masses are model features rather than distributed structural
    // elements. Include them only when requested by the load definition.
    if (consider_point_masses_) {
        apply_point_mass_inertial_contribution(model_data, bc, center_, center_acc_, omega_, alpha_);
    }
}

/**
 * Builds the diagnostic representation of the rigid-body inertia load.
 *
 * The result reports the target region, translational and angular kinematics,
 * and whether point-mass contributions are included.
 *
 * @return Human-readable load description.
 */
std::string InertialLoad::str() const {
    std::ostringstream os;

    // Report the region and every vector that determines the acceleration field,
    // followed by the separate point-mass participation flag.
    os << "INERTIAL: target=ELSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", center=[" << center_      (0) << ", " << center_      (1) << ", " << center_      (2) << "]"
       << ", a0=["     << center_acc_  (0) << ", " << center_acc_  (1) << ", " << center_acc_  (2) << "]"
       << ", omega=["  << omega_       (0) << ", " << omega_       (1) << ", " << omega_       (2) << "]"
       << ", alpha=["  << alpha_       (0) << ", " << alpha_       (1) << ", " << alpha_       (2) << "]"
       << ", consider_point_masses="
       << (consider_point_masses_ ? "true" : "false");

    return os.str();
}
} // namespace bc
} // namespace fem
