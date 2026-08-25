/**
 * @file load_inertial.cpp
 * @brief Implements equivalent forces for prescribed rigid-body acceleration.
 *
 * Distributed structural mass is handled through each element's density-scaled
 * vector-field integrator. Point elements are treated explicitly so their
 * concentrated translational mass and diagonal rotary inertia from
 * `PointMassSection` both contribute to rigid-body inertia loading.
 *
 * Regular compiled PointElements participate when their element id belongs to
 * the selected ELSET. Auxiliary point elements created by native
 * `POINTMASS, NSET=...` remain outside ELSET topology and are included only when
 * `consider_point_masses_` is enabled, preserving the existing FEMaster option.
 *
 * @see load_inertial.h
 * @see model::PointElement
 * @see PointMassSection
 * @author Finn Eggers
 * @date 25.08.2026
 */

#include "load_inertial.h"

#include "../core/logging.h"
#include "../model/element/element_structural.h"
#include "../model/element/point.h"
#include "../model/model_data.h"
#include "../section/section_point_mass.h"

#include <Eigen/Geometry>

#include <sstream>

namespace fem {
namespace bc {
namespace {

/**
 * Adds the concentrated rigid-body inertia force and moment of one PointElement.
 *
 * Translational acceleration contains reference-point, tangential and
 * centripetal terms. The section mass multiplies that acceleration directly at
 * the element node. Diagonal rotary inertia contributes Euler and gyroscopic
 * moments in the three nodal rotational directions.
 *
 * @param point Point element carrying a PointMassSection.
 * @param positions Global nodal position field.
 * @param bc Generalized nodal field receiving the result.
 * @param center Reference point of the rigid-body motion.
 * @param center_acc Translational acceleration at the reference point.
 * @param omega Angular velocity.
 * @param alpha Angular acceleration.
 */
void apply_point_element_inertial_contribution(const model::PointElement& point,
                                               const model::Field& positions,
                                               model::Field& bc,
                                               const Vec3& center,
                                               const Vec3& center_acc,
                                               const Vec3& omega,
                                               const Vec3& alpha) {
    if (!point._section) return;

    const auto* section = point._section->as<PointMassSection>();
    logging::error(section != nullptr,
        "InertialLoad: point element ", point.elem_id, " has a non-point-mass section");

    const ID node_id = point.nodes()[0];
    logging::error(node_id >= 0 && static_cast<Index>(node_id) < positions.rows,
        "InertialLoad: point-element node ", node_id,
        " is outside positions field with ", positions.rows, " rows");

    const Index node = static_cast<Index>(node_id);
    const Vec3 x     = positions.row_vec3(node);
    const Vec3 r     = x - center;
    const Vec3 a_rot = alpha.cross(r) + omega.cross(omega.cross(r));
    const Vec3 dF    = -section->mass_ * (center_acc + a_rot);

    bc(node, 0) += dF(0);
    bc(node, 1) += dF(1);
    bc(node, 2) += dF(2);

    if (bc.components >= 6) {
        const Vec3 Jalpha = section->rotary_inertia_.cwiseProduct(alpha);
        const Vec3 Jomega = section->rotary_inertia_.cwiseProduct(omega);
        const Vec3 dM     = -(Jalpha + omega.cross(Jomega));

        bc(node, 3) += dM(0);
        bc(node, 4) += dM(1);
        bc(node, 5) += dM(2);
    }
}

} // namespace

/**
 * Assembles equivalent inertia loads for distributed and concentrated mass.
 *
 * Regular elements in `region_` are evaluated first. A PointElement uses its
 * concentrated section directly; other structural elements integrate the
 * acceleration field with density scaling. When `consider_point_masses_` is set,
 * auxiliary post-compile PointElements created by native POINTMASS commands are
 * added independently because they do not belong to the dense ELSET namespace.
 *
 * @param model_data Model topology, fields and point-element storage.
 * @param bc Generalized nodal field receiving the contribution.
 * @param time Unused analysis time retained by the common interface.
 * @param ignore_amplitude Unused common-interface flag.
 */
void InertialLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude) {
    (void) time;
    (void) ignore_amplitude;

    logging::error(region_ != nullptr,
        "InertialLoad: region not set");
    logging::error(model_data.positions != nullptr,
        "InertialLoad: positions field not set in model data");

    const auto& positions = *model_data.positions;

    for (const ID el_id : *region_) {
        logging::error(el_id >= 0 && static_cast<std::size_t>(el_id) < model_data.elements.size(),
            "InertialLoad: element ", el_id, " is outside compiled element storage");

        auto& el_ptr = model_data.elements[static_cast<std::size_t>(el_id)];
        if (!el_ptr) continue;

        if (auto* point = el_ptr->as<model::PointElement>()) {
            apply_point_element_inertial_contribution(
                *point, positions, bc, center_, center_acc_, omega_, alpha_);
            continue;
        }

        auto* structural = el_ptr->as<model::StructuralElement>();
        if (!structural) continue;

        auto acceleration = [c = center_, a0 = center_acc_, w = omega_, al = alpha_](const Vec3& x) -> Vec3 {
            const Vec3 r = x - c;
            const Vec3 a_rot = al.cross(r) + w.cross(w.cross(r));
            return -(a0 + a_rot);
        };
        structural->integrate_vector_field(bc, true, acceleration);
    }

    if (consider_point_masses_) {
        for (const auto& element : model_data.point_elements) {
            if (!element) continue;
            const auto* point = element->as<model::PointElement>();
            logging::error(point != nullptr,
                "InertialLoad: auxiliary point-element storage contains a non-point element");
            apply_point_element_inertial_contribution(
                *point, positions, bc, center_, center_acc_, omega_, alpha_);
        }
    }
}

/**
 * Builds the diagnostic representation of the rigid-body inertia load.
 *
 * @return Human-readable load description.
 */
std::string InertialLoad::str() const {
    std::ostringstream os;

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
