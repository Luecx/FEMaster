/**
 * @file load_inertial.cpp
 * @brief Implements equivalent forces for prescribed rigid-body acceleration.
 *
 * Distributed structural mass is handled through each element's density-scaled
 * vector-field integrator. Point elements are treated explicitly so their
 * concentrated translational mass and diagonal rotary inertia from
 * `PointMassSection` both contribute to rigid-body inertia loading.
 *
 * Regular compiled point elements participate when their element identifier
 * belongs to the selected ELSET. Auxiliary point elements created by native
 * point-mass features remain outside compiled ELSET topology and are included
 * only when `consider_point_masses_` is enabled.
 *
 * @see load_inertial.h
 * @see model::PointElement
 * @see PointMassSection
 *
 * @author Finn Eggers
 * @date 25.08.2026
 */

#include "load_inertial.h"

#include "../../core/logging.h"
#include "../../model/element/element_structural.h"
#include "../../model/element/point.h"
#include "../../model/model_data.h"
#include "../../section/section_point_mass.h"

#include <Eigen/Geometry>
#include <sstream>

namespace fem::bc {
namespace {

/**
 * Adds the concentrated rigid-body inertia force and moment of one point element.
 *
 * Translational acceleration contains reference-point, tangential and
 * centripetal terms. The point-mass section multiplies that acceleration by its
 * concentrated mass and applies the opposite sign. Diagonal rotary inertia
 * contributes Euler and gyroscopic moments in the three nodal rotational
 * directions when the target field provides those components.
 *
 * @param point Point element carrying a `PointMassSection`.
 * @param positions Global nodal position field.
 * @param rhs Generalized nodal right-hand-side field receiving the contribution.
 * @param center Reference point of the prescribed rigid-body motion.
 * @param center_acc Translational acceleration at the reference point.
 * @param omega Angular velocity.
 * @param alpha Angular acceleration.
 */
void apply_point_element_inertial_contribution(const model::PointElement& point,
                                               const model::Field&        positions,
                                               model::Field&              rhs,
                                               const Vec3&                center,
                                               const Vec3&                center_acc,
                                               const Vec3&                omega,
                                               const Vec3&                alpha) {
    // A point element without a section carries no concentrated mass or rotary
    // inertia and therefore contributes nothing to the inertia load.
    if (!point._section) return;

    // Inertia loading of a point element requires the dedicated point-mass
    // section because generic point sections do not define mass properties.
    const auto* section = point._section->as<PointMassSection>();
    logging::error(section != nullptr,
        "InertialLoad: point element ", point.elem_id, " has a non-point-mass section");

    // Validate the single point-element node before accessing its physical
    // position in the compiled nodal field.
    const ID node_id = point.nodes()[0];
    logging::error(node_id >= 0 && static_cast<Index>(node_id) < positions.rows,
        "InertialLoad: point-element node ", node_id,
        " is outside positions field with ", positions.rows, " rows");

    // Evaluate rigid-body acceleration at the concentrated mass position. The
    // equivalent nodal inertia force acts opposite to the prescribed acceleration.
    const Index node  = static_cast<Index>(node_id);
    const Vec3  x     = positions.row_vec3(node);
    const Vec3  r     = x - center;
    const Vec3  a_rot = alpha.cross(r) + omega.cross(omega.cross(r));
    const Vec3  dF    = -section->mass_ * (center_acc + a_rot);

    rhs(node, 0) += dF(0);
    rhs(node, 1) += dF(1);
    rhs(node, 2) += dF(2);

    // Add rotary-inertia moments only when the target field carries structural
    // rotational components. Diagonal principal inertia is stored componentwise.
    if (rhs.components >= 6) {
        const Vec3 Jalpha = section->rotary_inertia_.cwiseProduct(alpha);
        const Vec3 Jomega = section->rotary_inertia_.cwiseProduct(omega);
        const Vec3 dM     = -(Jalpha + omega.cross(Jomega));

        rhs(node, 3) += dM(0);
        rhs(node, 4) += dM(1);
        rhs(node, 5) += dM(2);
    }
}

} // namespace

/**
 * Assembles equivalent inertia loads for distributed and concentrated mass.
 *
 * Regular elements in `region_` are evaluated first. A `PointElement` uses its
 * concentrated mass section directly; other structural elements integrate the
 * negative rigid-body acceleration field with density scaling. When
 * `consider_point_masses_` is set, auxiliary post-compile point elements are
 * added independently because they do not participate in normal ELSET topology.
 *
 * The current inertia definition is time-independent and therefore does not use
 * the common amplitude arguments.
 *
 * @param model_data Model topology, fields and auxiliary point-element storage.
 * @param rhs Generalized nodal right-hand-side field receiving the contribution.
 * @param time Unused analysis time retained by the common load interface.
 * @param ignore_amplitude Unused common-interface amplitude flag.
 */
void InertialLoad::apply(model::ModelData& model_data,
                         model::Field&     rhs,
                         Precision         time,
                         bool              ignore_amplitude) {
    // The rigid-body kinematics are stored directly on this load definition and
    // are not scaled by the common time/amplitude mechanism.
    (void) time;
    (void) ignore_amplitude;

    // Validate the selected mass region, nodal geometry and structural target
    // field before traversing compiled elements.
    logging::error(region_ != nullptr,
        "InertialLoad: region not set");
    logging::error(model_data.positions != nullptr,
        "InertialLoad: positions field not set in model data");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components >= 6,
        "InertialLoad: target field must be a NODE field with at least six components");

    const auto& positions = *model_data.positions;

    // Process all regular compiled elements selected by the element region.
    for (const ID element_id : *region_) {
        // Guard the dense compiled element storage before dereferencing the
        // region identifier.
        logging::error(element_id >= 0
                    && static_cast<std::size_t>(element_id) < model_data.elements.size(),
            "InertialLoad: element ", element_id, " is outside compiled element storage");

        auto& element = model_data.elements[static_cast<std::size_t>(element_id)];
        if (!element) continue;

        if (auto* point = element->as<model::PointElement>()) {
            // Concentrated point mass and rotary inertia are assembled directly
            // because no distributed element integration is required.
            apply_point_element_inertial_contribution(
                *point, positions, rhs, center_, center_acc_, omega_, alpha_);
            continue;
        }

        // Distributed inertia requires the structural vector-field integrator,
        // which supplies density, interpolation and the element measure.
        auto* structural = element->as<model::StructuralElement>();
        if (!structural) continue;

        // Construct the negative rigid-body acceleration field evaluated at
        // physical integration-point positions inside the structural element.
        auto acceleration = [c = center_, a0 = center_acc_, w = omega_, al = alpha_](const Vec3& x) -> Vec3 {
            const Vec3 r = x - c;
            return -(a0 + al.cross(r) + w.cross(w.cross(r)));
        };
        structural->integrate_vector_field(rhs, true, acceleration);
    }

    // Auxiliary point-mass features are stored outside the compiled element
    // region namespace. Include them explicitly only when requested by the load.
    if (consider_point_masses_) {
        for (const auto& element : model_data.point_elements) {
            if (!element) continue;

            const auto* point = element->as<model::PointElement>();
            logging::error(point != nullptr,
                "InertialLoad: auxiliary point-element storage contains a non-point element");

            apply_point_element_inertial_contribution(
                *point, positions, rhs, center_, center_acc_, omega_, alpha_);
        }
    }
}

/**
 * Builds the diagnostic representation of the rigid-body inertia load.
 *
 * @return Human-readable target, rigid-body kinematics and point-mass option.
 */
std::string InertialLoad::str() const {
    std::ostringstream os;

    // Report all kinematic parameters explicitly so the generated acceleration
    // field can be reconstructed from the model overview.
    os << "INERTIAL: target=ELSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", center=[" << center_(0) << ", " << center_(1) << ", " << center_(2) << "]"
       << ", a0=[" << center_acc_(0) << ", " << center_acc_(1) << ", " << center_acc_(2) << "]"
       << ", omega=[" << omega_(0) << ", " << omega_(1) << ", " << omega_(2) << "]"
       << ", alpha=[" << alpha_(0) << ", " << alpha_(1) << ", " << alpha_(2) << "]"
       << ", consider_point_masses=" << (consider_point_masses_ ? "true" : "false");

    return os.str();
}

} // namespace fem::bc
