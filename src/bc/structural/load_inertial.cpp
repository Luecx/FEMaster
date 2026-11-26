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
 * @see Condition
 * @see model::StructuralElement
 * @see model::PointElement
 * @see PointMassSection
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#include "load_inertial.h"

#include "../../core/logging.h"
#include "../../math/vec_util.h"
#include "../../model/element/element_structural.h"
#include "../../model/element/point.h"
#include "../../model/model_data.h"
#include "../../section/section_point_mass.h"

#include <Eigen/Geometry>

#include <sstream>

namespace fem::bc {

/**
 * Assembles equivalent inertia loads from prescribed rigid-body kinematics.
 *
 * A rigid-body point at x has acceleration
 *
 *     a(x) = a0 + alpha x (x-c) + omega x (omega x (x-c)).
 *
 * Distributed structural mass contributes -integral(rho N^T a dOmega) over the
 * reference volume, so the mass is conserved. Point masses contribute -m a(x),
 * and their optional diagonal rotary inertia contributes the Euler and
 * gyroscopic moments -J alpha - omega x (J omega).
 *
 * The calling analysis scales the complete inertia resultant. It does not
 * interpolate angular velocities, which enter some inertia terms quadratically.
 *
 * A dead load evaluates a(x) at the reference positions X. A follower load
 * evaluates it at the current positions x = X + u. Because the acceleration
 * field is affine in the position,
 *
 *     a(X + u) = a(X) + D u,   D = [alpha]x + omega omega^T - |omega|^2 I,
 *
 * the follower load of distributed mass is exactly
 *
 *     f = f_dead + K_L u,   K_L(i,j) = -s M_ij D,
 *
 * where M_ij is the density-scaled shape-function product operator of the
 * element field integration. It uses the same quadrature and nodal weighting
 * as f_dead, so the sum is the acceleration field integrated at the current
 * positions. The analytic load tangent df/du is therefore K_L itself. Point
 * masses evaluate a(x) directly at the current position with the tangent
 * -s m D. Rotary inertia moments do not depend on the position.
 *
 * Compiled point elements are selected through region_. Auxiliary point
 * elements are included separately if consider_point_masses_ is enabled.
 * A follower load without a requested load tangent is reported as an error.
 *
 * @param model_data Compiled structural topology, positions and point masses.
 * @param load Global nodal force and moment field modified in place.
 * @param load_tangent Optional load-tangent triplets receiving df/du.
 * @param system_dof_ids Global equation numbering of the load tangent.
 * @param scale Multiplier of the complete inertia contribution.
 */
void InertialLoad::apply(
    model::ModelData&        model_data,
    model::Field&            load,
    TripletList*             load_tangent,
    constraint::Equations&,
    const model::SystemDofs& system_dof_ids,
    TripletList&,
    Precision                scale
) {
    // Validate the definition and the configuration required by the flag
    logging::error(region_ != nullptr,
        "INERTIAL: target element region is not initialized");
    logging::error(model_data.positions != nullptr && model_data.positions_reference != nullptr,
        "INERTIAL: positions fields are not initialized");
    logging::error(!follower_ || load_tangent != nullptr,
        "INERTIAL: follower load requires a load tangent and is only supported in nonlinear analyses");

    if (scale == Precision(0)) return;

    // Dead loads use the reference positions, follower loads the current ones
    const auto& reference = *model_data.positions_reference;
    const auto& positions = follower_ ? *model_data.positions : reference;

    // Evaluate complete centripetal and Euler acceleration before scaling.
    // Scaling omega itself would produce an incorrect quadratic contribution.
    const auto inertia_acceleration = [&](const Vec3& position) -> Vec3 {
        const Vec3 r = position - center_;
        const Vec3 acceleration =
            center_acc_ + alpha_.cross(r) + omega_.cross(omega_.cross(r));
        return -scale * acceleration;
    };

    // Derivative of the scaled inertia acceleration -s a(x) with respect to x
    const Mat3 acceleration_gradient = -scale * (
        math::skew(alpha_)
      + omega_ * omega_.transpose()
      - omega_.squaredNorm() * Mat3::Identity());

    // Apply the concentrated inertia load to a compiled or auxiliary point
    // element. Both storage paths use the same force and moment calculation.
    const auto apply_point_mass = [&](const model::PointElement& point) {
        if (!point._section) return;

        const ID    node_id = point.nodes()[0];
        const auto* section = point._section->as<PointMassSection>();

        logging::error(section != nullptr,
            "INERTIAL: point element ", point.label(), " has a non-point-mass section");
        logging::error(node_id >= 0 && static_cast<Index>(node_id) < positions.rows,
            "INERTIAL: point-element node ", node_id, " is outside positions field with ", positions.rows, " rows");

        const Index node     = static_cast<Index>(node_id);
        const Vec3  position = positions.row_vec3(node);

        // F = -m a(x), including the effective step multiplier
        const Vec3 force = section->mass_ * inertia_acceleration(position);

        for (Dim i = 0; i < 3; ++i) {
            load(node, i) += force[i];
        }

        // The concentrated force follows the current position of the node
        if (follower_) {
            const Mat3 stiffness = section->mass_ * acceleration_gradient;
            add_load_tangent_block(*load_tangent, system_dof_ids, node_id, 0, node_id, stiffness);
        }

        if (load.components < 6) return;

        // Scale the complete Euler and gyroscopic moments at this point mass
        const Vec3 J_omega = section->rotary_inertia_.cwiseProduct(omega_);
        const Vec3 moment = -scale * (
            section->rotary_inertia_.cwiseProduct(alpha_) + omega_.cross(J_omega));

        for (Dim i = 0; i < 3; ++i) {
            load(node, i + 3) += moment[i];
        }
    };

    // Add the follower part K_L u of distributed mass. The shape-function
    // products use the same reference quadrature and nodal weighting as the
    // dead field integration, so both parts form one consistent integral.
    const auto apply_follower_mass = [&](model::StructuralElement& element) {
        const Index         n_nodes  = static_cast<Index>(element.n_nodes());
        const ID*           nodes    = element.nodes();
        const DynamicMatrix products = element.integrate_shape_products(true);

        logging::error(products.rows() == n_nodes && products.cols() == n_nodes,
            "INERTIAL: element ", element.label(), " returned an invalid shape-product operator");

        for (Index i = 0; i < n_nodes; ++i) {
            for (Index j = 0; j < n_nodes; ++j) {
                if (products(i, j) == Precision(0)) continue;

                // Exact follower increment and its derivative
                const Mat3 stiffness    = products(i, j) * acceleration_gradient;
                const Vec3 displacement = positions.row_vec3(nodes[j]) - reference.row_vec3(nodes[j]);
                const Vec3 force        = stiffness * displacement;

                for (Dim d = 0; d < 3; ++d) {
                    load(nodes[i], d) += force[d];
                }

                add_load_tangent_block(*load_tangent, system_dof_ids, nodes[i], 0, nodes[j], stiffness);
            }
        }
    };

    // Distributed mass is integrated with density scaling over the reference
    // volume. Compiled PointElements use the concentrated-mass formula instead.
    for (ID element_id : *region_) {
        logging::error(element_id >= 0 && static_cast<std::size_t>(element_id) < model_data.elements.size(),
            "INERTIAL: element ", element_id, " is outside compiled element storage");

        const auto& element = model_data.elements[static_cast<std::size_t>(element_id)];
        if (!element) continue;

        if (const auto* point = element->as<model::PointElement>()) {
            apply_point_mass(*point);
        } else if (auto* structural = element->as<model::StructuralElement>()) {
            structural->integrate_vector_field(load, true, inertia_acceleration);

            // A purely translational acceleration such as gravity is uniform in
            // space, so its follower part K_L u vanishes identically
            if (follower_ && !acceleration_gradient.isZero()) apply_follower_mass(*structural);
        }
    }

    if (!consider_point_masses_) return;

    // Native POINTMASS definitions created after compilation do not belong
    // to an ELSET and must be assembled in a separate pass.
    for (const auto& element : model_data.point_elements) {
        if (!element) continue;

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
