/**
 * @file frt_shell_integrate.cpp
 * @brief Implements shell mass, volume and distributed-field integration.
 *
 * All integrations use the actual isoparametric reference midsurface Jacobian.
 * Translational and rotational consistent mass contributions are integrated at
 * every quadrature point. The rotational inertia projector follows the
 * pointwise normalized reference director rather than one constant planar
 * element normal.
 *
 * Follower-pressure terms are intentionally not part of this file. Distributed
 * scalar, vector and tensor fields use the reference volume measure and are
 * evaluated at the current midsurface position.
 *
 * @see FRTShell
 *
 * @author Finn Eggers
 * @date 20.07.2026
 */

#include "frt_shell.h"

#include "../../core/logging.h"

namespace fem::model {

template<Index N>
Precision FRTShell<N>::volume() {
    return this->get_section()->thickness_ * reference_data().area;
}

/**
 * Integrates the scalar shape-function product matrix over the reference
 * midsurface.
 */
template<Index N>
StaticMatrix<N, N> FRTShell<N>::integrate_NNt() const {
    StaticMatrix<N, N> result;
    result.setZero();

    for (const ReferencePoint& point : reference_data().ip_points) {
        result.noalias() += point.w
                          * point.detJ
                          * point.shape
                          * point.shape.transpose();
    }

    return result;
}

/**
 * Assembles the consistent translational and rotational shell mass matrix.
 *
 * Translational inertia uses `rho h N_i N_j`. Rotational inertia uses
 * `rho h^3 / 12` projected into the local tangent plane of the interpolated
 * reference director. A very small drilling inertia prevents a singular mass
 * matrix for the regularized drilling degree of freedom.
 */
template<Index N>
MapMatrix FRTShell<N>::mass(Precision* buffer) {
    const Precision rho = this->get_density(false);

    Mat6N mass_matrix;
    mass_matrix.setZero();

    if (rho > Precision(0)) {
        const Precision h           = this->get_section()->thickness_;
        const Precision mass_area   = rho * h;
        const Precision inertia_rot = rho * h * h * h / Precision(12);

        // Integrate the pointwise rotational projector on the curved reference
        // surface instead of using one constant element normal
        for (const ReferencePoint& point : reference_data().ip_points) {
            Vec3 director = point.D;

            if (director.squaredNorm() <= Precision(1e-24)) {
                director = point.e3;
            } else {
                director.normalize();
            }

            const Mat3 tangent_projector = Mat3::Identity() - director * director.transpose();
            const Mat3 drill_projector   = director * director.transpose();
            const Precision weight       = point.w * point.detJ;

            for (Index i = 0; i < num_nodes; ++i) {
                for (Index j = 0; j < num_nodes; ++j) {
                    const Precision Nij = weight * point.shape(i) * point.shape(j);

                    mass_matrix.template block<3, 3>(
                        dofs_per_node * i,
                        dofs_per_node * j
                    ) += mass_area * Nij * Mat3::Identity();

                    mass_matrix.template block<3, 3>(
                        dofs_per_node * i + 3,
                        dofs_per_node * j + 3
                    ) += inertia_rot * Nij * tangent_projector;

                    // Keep the drilling inertia several orders of magnitude
                    // below the physical rotational inertia
                    if (i == j) {
                        mass_matrix.template block<3, 3>(
                            dofs_per_node * i + 3,
                            dofs_per_node * i + 3
                        ) += Precision(1e-6) * mass_area * Nij * drill_projector;
                    }
                }
            }
        }
    }

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = mass_matrix;
    return mapped;
}

template<Index N>
Vec3 FRTShell<N>::global_point_current(Precision r, Precision s) const {
    const VecN  shape     = shape_function(r, s);
    const MatN6 positions = node_coords_current_6();

    Vec3 point = Vec3::Zero();

    for (Index node = 0; node < num_nodes; ++node) {
        point += shape(node) * positions.template block<1, 3>(node, 0).transpose();
    }

    return point;
}

/**
 * Integrates a scalar field through the shell volume.
 *
 * The field is evaluated at the current midsurface position. The integration
 * measure is `h dA0`, optionally multiplied by the material density.
 */
template<Index N>
Precision FRTShell<N>::integrate_scalar_field(bool               scale_by_density,
                                              const ScalarField& field) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

    Precision result = Precision(0);

    for (const ReferencePoint& point : reference_data().ip_points) {
        result += field(global_point_current(point.r, point.s))
                * rho
                * h
                * point.w
                * point.detJ;
    }

    return result;
}

/**
 * Integrates a vector field through the shell volume.
 */
template<Index N>
Vec3 FRTShell<N>::integrate_vector_field(bool            scale_by_density,
                                         const VecField& field) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

    Vec3 result = Vec3::Zero();

    for (const ReferencePoint& point : reference_data().ip_points) {
        result += field(global_point_current(point.r, point.s))
                * (rho * h * point.w * point.detJ);
    }

    return result;
}

/**
 * Integrates a distributed vector field and scatters its consistent nodal load.
 */
template<Index N>
void FRTShell<N>::integrate_vector_field(Field&          node_loads,
                                         bool            scale_by_density,
                                         const VecField& field) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

    for (const ReferencePoint& point : reference_data().ip_points) {
        const Vec3 force = field(global_point_current(point.r, point.s))
                         * (rho * h * point.w * point.detJ);

        for (Index node = 0; node < num_nodes; ++node) {
            const Index node_id = static_cast<Index>(this->node_ids[node]);

            for (Index component = 0; component < 3; ++component) {
                node_loads(node_id, component) += point.shape(node) * force(component);
            }
        }
    }
}

/**
 * Integrates a second-order tensor field through the shell volume.
 */
template<Index N>
Mat3 FRTShell<N>::integrate_tensor_field(bool            scale_by_density,
                                         const TenField& field) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

    Mat3 result = Mat3::Zero();

    for (const ReferencePoint& point : reference_data().ip_points) {
        result += field(global_point_current(point.r, point.s))
                * (rho * h * point.w * point.detJ);
    }

    return result;
}

template<Index N>
void FRTShell<N>::compute_compliance(Field& displacement, Field& result) {
    Precision buffer[num_dofs * num_dofs];
    MapMatrix K = stiffness(buffer);
    const Vec6N u = element_displacement_vector(displacement);

    result(static_cast<Index>(this->elem_id), 0) = u.dot(K * u);
}

template<Index N>
void FRTShell<N>::compute_compliance_angle_derivative(Field& displacement,
                                                      Field& result) {
    (void) displacement;
    (void) result;
}

template<Index N>
bool FRTShell<N>::compute_shear_flow(Field&       shear_flow,
                                     const Field& displacement,
                                     int          offset) {
    (void) shear_flow;
    (void) displacement;
    (void) offset;
    return false;
}

template<Index N>
bool FRTShell<N>::compute_beam_section_forces(Field&       section_forces,
                                              const Field& displacement,
                                              int          offset) {
    (void) section_forces;
    (void) displacement;
    (void) offset;
    return false;
}

} // namespace fem::model
