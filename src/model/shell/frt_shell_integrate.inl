/**
 * @file frt_shell_integrate.inl
 * @brief Implements shell mass, volume and distributed-field integration.
 *
 * All integrations use the actual isoparametric reference midsurface Jacobian.
 * Translational and rotational consistent mass contributions are integrated at
 * every quadrature point. The rotational inertia projector follows the
 * pointwise normalized reference director instead of one planar element normal.
 *
 * Distributed fields are evaluated at the current midsurface position while
 * the integration measure remains the Total-Lagrangian reference volume.
 *
 * @see FRTShell
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#include "frt_shell.h"

namespace fem::model {

/**
 * Returns the reference volume of the shell element.
 *
 * The shell thickness is assumed constant over the element and multiplies the
 * integrated curved reference midsurface area.
 *
 * @return Reference shell volume.
 */
template<Index N>
Precision FRTShell<N>::volume() {
    return this->get_section()->thickness_ * reference_data().area;
}

/**
 * Assembles the consistent translational and rotational shell mass matrix.
 *
 * Translational inertia uses `rho h N_i N_j`. Physical rotary inertia uses
 * `rho h^3/12` in the tangent plane normal to the interpolated reference
 * director. A dimensionless fraction of the same physical rotary inertia is
 * assigned to the drilling direction so the regularized rotational mass block
 * remains nonsingular without introducing the previous dimensionally incorrect
 * `rho h` drilling term.
 *
 * @param buffer Caller-provided dense element matrix storage.
 * @return Mapped consistent shell mass matrix.
 */
template<Index N>
MapMatrix FRTShell<N>::mass(Precision* buffer) {
    const Precision rho = this->get_density(false);

    Mat6N mass_matrix = Mat6N::Zero();

    if (rho > Precision(0)) {
        const Precision h           = this->get_section()->thickness_;
        const Precision mass_area   = rho * h;
        const Precision inertia_rot = rho * h * h * h / Precision(12);

        for (const ReferencePoint& point : reference_data().ip_points) {
            Vec3 director = point.D;

            // Use the interpolated reference director where possible and fall
            // back to the exact surface normal for a degenerate interpolation
            if (director.squaredNorm() <= Precision(1e-24)) {
                director = point.e3;
            } else {
                director.normalize();
            }

            const Mat3 drill_projector   = director * director.transpose();
            const Mat3 tangent_projector = Mat3::Identity() - drill_projector;
            const Mat3 rotation_metric   = tangent_projector
                                         + drill_inertia_scale * drill_projector;
            const Precision weight       = point.w * point.detJ;

            // Assemble both translational and rotational consistent N_i N_j
            // blocks for every nodal pair
            for (Index i = 0; i < num_nodes; ++i) {
                for (Index j = 0; j < num_nodes; ++j) {
                    const Precision Nij = weight * point.shape(i) * point.shape(j);

                    mass_matrix.template block<3, 3>(
                        dofs_per_node * i,
                        dofs_per_node * j
                    ).diagonal().array() += mass_area * Nij;

                    mass_matrix.template block<3, 3>(
                        dofs_per_node * i + 3,
                        dofs_per_node * j + 3
                    ) += inertia_rot * Nij * rotation_metric;
                }
            }
        }
    }

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = mass_matrix;
    return mapped;
}

/**
 * Interpolates the current physical midsurface position at one natural point.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @return Current global midsurface position.
 */
template<Index N>
Vec3 FRTShell<N>::global_point_current(Precision r, Precision s) const {
    const VecN  shape     = shape_function(r, s);
    const MatN6 positions = node_coords_current_6();

    Vec3 point = Vec3::Zero();

    // Interpolate only the three translational position components
    for (Index node = 0; node < num_nodes; ++node) {
        point += shape(node)
               * positions.template block<1, 3>(node, 0).transpose();
    }

    return point;
}

/**
 * Integrates a scalar field through the shell reference volume.
 *
 * Current nodal positions are gathered once and reused at every integration
 * point. The field is evaluated at the current midsurface position, while the
 * integration measure is `h dA0` and may optionally include material density.
 *
 * @param scale_by_density Multiply the integrand by the current material density.
 * @param field Spatial scalar field evaluated in global coordinates.
 * @return Integrated scalar quantity.
 */
template<Index N>
Precision FRTShell<N>::integrate_scalar_field(
    bool               scale_by_density,
    const ScalarField& field
) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);
    const MatN6 positions = node_coords_current_6();

    Precision result = Precision(0);

    for (const ReferencePoint& point : reference_data().ip_points) {
        Vec3 position = Vec3::Zero();

        // Interpolate the current physical evaluation position from the nodal
        // translations already gathered before the quadrature loop
        for (Index node = 0; node < num_nodes; ++node) {
            position += point.shape(node)
                      * positions.template block<1, 3>(node, 0).transpose();
        }

        const Precision weighted_volume =
            rho * h * point.w * point.detJ;
        result += field(position) * weighted_volume;
    }

    return result;
}

/**
 * Integrates a vector field through the shell reference volume.
 *
 * @param scale_by_density Multiply the integrand by the current material density.
 * @param field Spatial vector field evaluated in global coordinates.
 * @return Integrated global vector.
 */
template<Index N>
Vec3 FRTShell<N>::integrate_vector_field(
    bool            scale_by_density,
    const VecField& field
) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);
    const MatN6 positions = node_coords_current_6();

    Vec3 result = Vec3::Zero();

    for (const ReferencePoint& point : reference_data().ip_points) {
        Vec3 position = Vec3::Zero();

        // Reuse the gathered nodal positions for the complete quadrature loop
        for (Index node = 0; node < num_nodes; ++node) {
            position += point.shape(node)
                      * positions.template block<1, 3>(node, 0).transpose();
        }

        const Precision weighted_volume =
            rho * h * point.w * point.detJ;
        result += field(position) * weighted_volume;
    }

    return result;
}

/**
 * Integrates a distributed vector field and scatters its consistent nodal load.
 *
 * The current physical field value is multiplied by `N_i h dA0` and optional
 * density before being added to every translational nodal degree of freedom.
 *
 * @param node_loads Global nodal load field to increment.
 * @param scale_by_density Multiply the prescribed field by material density.
 * @param field Spatial vector field evaluated in global coordinates.
 */
template<Index N>
void FRTShell<N>::integrate_vector_field(
    Field&          node_loads,
    bool            scale_by_density,
    const VecField& field
) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);
    const MatN6 positions = node_coords_current_6();

    for (const ReferencePoint& point : reference_data().ip_points) {
        Vec3 position = Vec3::Zero();

        // Evaluate the current global point associated with this quadrature
        // coordinate
        for (Index node = 0; node < num_nodes; ++node) {
            position += point.shape(node)
                      * positions.template block<1, 3>(node, 0).transpose();
        }

        const Vec3 force = field(position)
                         * (rho * h * point.w * point.detJ);

        // Distribute the integrated field consistently through the same shape
        // functions used by the shell geometry
        for (Index node = 0; node < num_nodes; ++node) {
            const Index node_id = static_cast<Index>(this->node_ids[node]);

            for (Index component = 0; component < 3; ++component) {
                node_loads(node_id, component) +=
                    point.shape(node) * force(component);
            }
        }
    }
}

/**
 * Integrates a second-order tensor field through the shell reference volume.
 *
 * @param scale_by_density Multiply the integrand by the current material density.
 * @param field Spatial tensor field evaluated in global coordinates.
 * @return Integrated global second-order tensor.
 */
template<Index N>
Mat3 FRTShell<N>::integrate_tensor_field(
    bool            scale_by_density,
    const TenField& field
) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);
    const MatN6 positions = node_coords_current_6();

    Mat3 result = Mat3::Zero();

    for (const ReferencePoint& point : reference_data().ip_points) {
        Vec3 position = Vec3::Zero();

        // Interpolate the physical field coordinate without regathering nodal
        // data inside the quadrature loop
        for (Index node = 0; node < num_nodes; ++node) {
            position += point.shape(node)
                      * positions.template block<1, 3>(node, 0).transpose();
        }

        const Precision weighted_volume =
            rho * h * point.w * point.detJ;
        result += field(position) * weighted_volume;
    }

    return result;
}

/**
 * Computes the element compliance contribution `u^T K u`.
 *
 * @param displacement Global nodal displacement field.
 * @param result Element result field receiving the scalar compliance.
 */
template<Index N>
void FRTShell<N>::compute_compliance(Field& displacement, Field& result) {
    Precision buffer[num_dofs * num_dofs];
    MapMatrix K = stiffness(buffer);
    const Vec6N u = element_displacement_vector(displacement);

    result(static_cast<Index>(this->elem_id), 0) = u.dot(K * u);
}

/**
 * Leaves the unsupported compliance derivative with respect to material angle
 * unchanged.
 *
 * @param displacement Unused global displacement field.
 * @param result Unused result field.
 */
template<Index N>
void FRTShell<N>::compute_compliance_angle_derivative(
    Field& displacement,
    Field& result
) {
    (void) displacement;
    (void) result;
}

/**
 * Reports that generic beam-style shear-flow recovery is unavailable for shell
 * elements.
 *
 * @param shear_flow Unused output field.
 * @param displacement Unused displacement field.
 * @param offset Unused output offset.
 * @return Always `false`.
 */
template<Index N>
bool FRTShell<N>::compute_shear_flow(
    Field&       shear_flow,
    const Field& displacement,
    int          offset
) {
    (void) shear_flow;
    (void) displacement;
    (void) offset;
    return false;
}

/**
 * Reports that generic beam-section force recovery is unavailable for shell
 * elements.
 *
 * @param section_forces Unused output field.
 * @param displacement Unused displacement field.
 * @param offset Unused output offset.
 * @return Always `false`.
 */
template<Index N>
bool FRTShell<N>::compute_beam_section_forces(
    Field&       section_forces,
    const Field& displacement,
    int          offset
) {
    (void) section_forces;
    (void) displacement;
    (void) offset;
    return false;
}

} // namespace fem::model
