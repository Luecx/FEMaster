/**
 * @file load_inertial.cpp
 * @brief Implements equivalent loads for prescribed rigid-body acceleration.
 *
 * For a point at global position `x`, the represented rigid-body kinematics use
 * the lever arm
 *
 *     r = x - c
 *
 * from the prescribed reference center `c` and the acceleration field
 *
 *     a(x) = a0 + alpha x r + omega x (omega x r).
 *
 * The second and third terms are the tangential and centripetal accelerations.
 * No Coriolis term appears because this load describes material points fixed in
 * the prescribed rigid-body motion rather than particles with relative velocity
 * in the rotating frame.
 *
 * Distributed structural mass contributes
 *
 *     f_e = -integral_Omega_e rho N^T a(x) dOmega.
 *
 * Point masses are handled explicitly so both concentrated translational mass
 * and diagonal rotary inertia can contribute to the generalized nodal load.
 *
 * @see InertialLoad
 * @see Neumann
 * @see model::StructuralElement
 * @see model::PointElement
 * @see PointMassSection
 *
 * @author Finn Eggers
 * @date 17.09.2026
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

/**
 * Assembles equivalent inertia loads for distributed and concentrated mass.
 *
 * The optional amplitude is evaluated once and multiplies the complete inertia
 * contribution. For distributed structural mass, the element integrator is
 * called with density scaling enabled so its physical measure becomes
 * `rho dOmega`. The supplied vector field is the negative rigid-body
 * acceleration,
 *
 *     b_inertia(x) = -a(t) [a0 + alpha x r + omega x (omega x r)],
 *
 * which yields the consistent element force
 *
 *     f_e = integral_Omega_e rho N^T b_inertia dOmega.
 *
 * A `PointElement` cannot be represented by a volume integral. Its translational
 * contribution is therefore assembled directly as
 *
 *     F = -a(t) m a(x).
 *
 * If the nodal RHS exposes rotational DOFs and the point-mass section stores
 * diagonal rotary inertia `J`, the Euler and gyroscopic moments are
 *
 *     M = -a(t) [J alpha + omega x (J omega)].
 *
 * Regular compiled point elements are selected through `region_`. Auxiliary
 * point elements created by native post-compile `POINTMASS` definitions are
 * traversed separately only when `consider_point_masses_` is enabled.
 *
 * @param model_data Compiled structural topology, nodal geometry and auxiliary
 *                   point-element storage.
 * @param rhs Generalized nodal RHS field modified in place.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Apply the nominal inertia load with unit amplitude
 *                         when true.
 */
void InertialLoad::apply(model::ModelData& model_data, model::Field& rhs, Precision time, bool ignore_amplitude) {
    // Validate the selected mass region and geometric field before constructing
    // the rigid-body acceleration
    logging::error(region_ != nullptr,
        "INERTIAL: target element region is not initialized");
    logging::error(model_data.positions != nullptr,
        "INERTIAL: positions field is not initialized");

    const auto&     positions = *model_data.positions;
    const Precision scale     = amplitude_ && !ignore_amplitude ? amplitude_->evaluate(time) : Precision(1);

    // Construct the signed acceleration field used by both distributed and
    // concentrated translational mass contributions. The returned vector already
    // includes the minus sign of the d'Alembert inertia force and the amplitude.
    auto inertia_acceleration = [&](const Vec3& position) -> Vec3 {
        const Vec3 radius                  = position - center_;
        const Vec3 tangential_acceleration = alpha_.cross(radius);
        const Vec3 centripetal_acceleration = omega_.cross(omega_.cross(radius));
        const Vec3 rigid_acceleration      = center_acc_ + tangential_acceleration + centripetal_acceleration;

        return -scale * rigid_acceleration;
    };

    // Assemble one concentrated PointElement. A local lambda keeps this
    // point-mass-specific operation beside the only algorithm that uses it while
    // allowing the regular and auxiliary point-element paths to share it.
    auto apply_point_mass = [&](const model::PointElement& point) {
        if (!point._section) {
            return;
        }

        const auto* section = point._section->as<PointMassSection>();
        logging::error(section != nullptr,
            "INERTIAL: point element ", point.elem_id, " has a non-point-mass section");

        const ID node_id = point.nodes()[0];
        logging::error(node_id >= 0 && static_cast<Index>(node_id) < positions.rows,
            "INERTIAL: point-element node ", node_id,
            " is outside positions field with ", positions.rows, " rows");

        const Index node     = static_cast<Index>(node_id);
        const Vec3  position = positions.row_vec3(node);

        // Translational concentrated inertia follows directly from F = m b,
        // where b is the signed inertia acceleration defined above
        const Vec3 force = section->mass_ * inertia_acceleration(position);
        rhs(node, 0) += force(0);
        rhs(node, 1) += force(1);
        rhs(node, 2) += force(2);

        if (rhs.components < 6) {
            return;
        }

        // The stored rotary inertia is diagonal in the nodal rotational basis.
        // Form J alpha and J omega componentwise, then add the Euler and
        // gyroscopic terms with the same d'Alembert sign and amplitude scale.
        const Vec3 inertia_alpha = section->rotary_inertia_.cwiseProduct(alpha_);
        const Vec3 inertia_omega = section->rotary_inertia_.cwiseProduct(omega_);
        const Vec3 moment        = -scale * (inertia_alpha + omega_.cross(inertia_omega));

        rhs(node, 3) += moment(0);
        rhs(node, 4) += moment(1);
        rhs(node, 5) += moment(2);
    };

    // Assemble the mass carried by regular compiled elements selected through
    // the element region
    for (ID element_id : *region_) {
        logging::error(element_id >= 0 && static_cast<std::size_t>(element_id) < model_data.elements.size(),
            "INERTIAL: element ", element_id, " is outside compiled element storage");

        auto& element = model_data.elements[static_cast<std::size_t>(element_id)];
        if (!element) {
            continue;
        }

        // A compiled PointElement represents concentrated rather than distributed
        // mass and must therefore use the direct point-mass equations
        if (auto* point = element->as<model::PointElement>()) {
            apply_point_mass(*point);
            continue;
        }

        // Distributed structural mass evaluates integral rho N^T b dOmega.
        // Density scaling is enabled explicitly because b is an acceleration.
        auto* structural = element->as<model::StructuralElement>();
        if (!structural) {
            continue;
        }

        structural->integrate_vector_field(rhs, true, inertia_acceleration);
    }

    if (!consider_point_masses_) {
        return;
    }

    // Native POINTMASS definitions created after compilation live outside the
    // dense element/ELSET namespace and are therefore added in a separate pass
    for (const auto& element : model_data.point_elements) {
        if (!element) {
            continue;
        }

        const auto* point = element->as<model::PointElement>();
        logging::error(point != nullptr,
            "INERTIAL: auxiliary point-element storage contains a non-point element");

        apply_point_mass(*point);
    }
}

/**
 * Builds a compact diagnostic representation of the rigid-body inertia load.
 *
 * The output reports the selected element region, reference center,
 * translational acceleration, angular velocity, angular acceleration and the
 * auxiliary-point-mass flag. Values are shown before amplitude scaling.
 *
 * @return Human-readable inertia-load description.
 */
std::string InertialLoad::str() const {
    std::ostringstream os;

    os << "INERTIAL: target=ELSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", center=[" << center_    (0) << ", " << center_    (1) << ", " << center_    (2) << "]"
       << ", a0=["     << center_acc_(0) << ", " << center_acc_(1) << ", " << center_acc_(2) << "]"
       << ", omega=["  << omega_     (0) << ", " << omega_     (1) << ", " << omega_     (2) << "]"
       << ", alpha=["  << alpha_     (0) << ", " << alpha_     (1) << ", " << alpha_     (2) << "]"
       << ", consider_point_masses=" << (consider_point_masses_ ? "true" : "false");

    if (amplitude_) {
        os << ", amplitude=" << amplitude_->name;
    }

    return os.str();
}

} // namespace fem::bc
