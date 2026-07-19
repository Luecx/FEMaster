/**
 * @file element_solid.ipp
 * @brief Implementation of the SolidElement class template. This file contains
 * the definitions for methods declared in the SolidElement class, including
 * the computation of the strain-displacement matrix, Jacobian, stiffness and
 * mass matrices, and other element-level calculations.
 *
 * @date Created on 12.06.2023
 * @author Finn Eggers
 */

#pragma once

#include "../../cos/rectangular_system.h"

namespace fem::model {

template<Index N>
void SolidElement<N>::compute_stress_strain(Field* strain,
                                            Field* stress,
                                            const Field& displacement,
                                            const RowMatrix& rst,
                                            int offset,
                                            bool use_green_lagrange_nl) {
    logging::error(strain != nullptr || stress != nullptr,
                   "SolidElement: compute_stress_strain requires at least one output field");
    logging::error(rst.cols() >= 3,
                   "SolidElement: stress/strain evaluation coordinates require at least 3 columns");

    auto reference_coords = this->node_coords_reference();
    auto current_coords = this->node_coords_current();
    auto local_displacement = this->nodal_data<3>(displacement);
    auto local_disp_mat = StaticMatrix<3, N>(local_displacement.transpose());
    auto local_displacement_vec = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    for (Eigen::Index n = 0; n < rst.rows(); ++n) {
        const Precision r = rst(n, 0);
        const Precision s = rst(n, 1);
        const Precision t = rst(n, 2);
        const Index row = static_cast<Index>(offset + n);

        if (!use_green_lagrange_nl) {
            Precision det;
            StaticMatrix<n_strain, D * N> B =
                this->strain_displacements(reference_coords, r, s, t, det, false);
            if (det <= Precision(0) || std::isnan(det) || std::isinf(det)) {
                continue;
            }

            const Vec6                   global_strain_voigt = B * local_displacement_vec;
            const VolumeStrainLinearized global_strain(global_strain_voigt);
            VolumeStressCauchy           global_stress;
            Mat6                         global_tangent;
            evaluate_material(r, s, t, global_strain, global_stress, global_tangent);

            for (Dim j = 0; j < n_strain; ++j) {
                if (strain) (*strain)(row, j) = global_strain.voigt()(j);
                if (stress) (*stress)(row, j) = global_stress.voigt()(j);
            }
            continue;
        }

        const Mat3 F = this->deformation_gradient(reference_coords, current_coords, r, s, t);

        const VolumeStrainGreenLagrange green_lagrange =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);
        VolumeStressPK2 second_pk;
        Mat6            tangent;
        evaluate_material(r, s, t, green_lagrange, second_pk, tangent);

        const VolumeStressCauchy cauchy = second_pk.to_cauchy(F);

        for (Dim j = 0; j < n_strain; ++j) {
            if (strain) (*strain)(row, j) = green_lagrange.voigt()(j);
            if (stress) (*stress)(row, j) = cauchy.voigt()(j);
        }
    }
}

template<Index N>
void SolidElement<N>::compute_stress_state(Field& stress_state,
                                           const Field& displacement,
                                           int offset,
                                           bool use_green_lagrange_nl) {
    const RowMatrix rst = stress_strain_ip_rst();
    if (rst.rows() == 0) {
        return;
    }

    const auto reference_coords       = this->node_coords_reference();
    const auto current_coords         = this->node_coords_current();
    const auto local_displacement     = this->nodal_data<3>(displacement);
    const auto local_disp_mat         = StaticMatrix<3, N>(local_displacement.transpose());
    const auto local_displacement_vec =
        Eigen::Map<const StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    for (Eigen::Index n = 0; n < rst.rows(); ++n) {
        const Precision r   = rst(n, 0);
        const Precision s   = rst(n, 1);
        const Precision t   = rst(n, 2);
        const Index     row = static_cast<Index>(offset + n);

        if (!use_green_lagrange_nl) {
            Precision det = Precision(0);
            const StaticMatrix<n_strain, D * N> B =
                this->strain_displacements(reference_coords, r, s, t, det, false);

            if (det <= Precision(0) || !std::isfinite(det)) {
                continue;
            }

            const Vec6                   strain_voigt = B * local_displacement_vec;
            const VolumeStrainLinearized strain(strain_voigt);
            VolumeStressCauchy           cauchy;
            Mat6                         tangent;
            evaluate_material(r, s, t, strain, cauchy, tangent);

            for (Dim component = 0; component < n_strain; ++component) {
                stress_state(row, component) = cauchy.voigt()(component);
            }
            continue;
        }

        const Mat3 F = this->deformation_gradient(reference_coords, current_coords, r, s, t);
        const VolumeStrainGreenLagrange green_lagrange =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);
        VolumeStressPK2 second_pk;
        Mat6            tangent;
        evaluate_material(r, s, t, green_lagrange, second_pk, tangent);

        // Total-Lagrange internal forces and geometric stiffness require PK2
        for (Dim component = 0; component < n_strain; ++component) {
            stress_state(row, component) = second_pk.voigt()(component);
        }
    }
}

template<Index N>
void SolidElement<N>::compute_internal_force_nonlinear(Field& node_forces,
                                                       const Field& ip_stress) {
    auto reference_coords = this->node_coords_reference();
    auto current_coords   = this->node_coords_current();
    auto scheme           = this->integration_scheme();

    for (Index n = 0; n < scheme.count(); n++) {
        const Precision r = scheme.get_point(n).r;
        const Precision s = scheme.get_point(n).s;
        const Precision t = scheme.get_point(n).t;
        const Precision w = scheme.get_point(n).w;
        const Index row = this->ip_index(n);

        Precision det0;
        // derivative of shape function w.r.t reference coordinates XYZ (N x 3)
        const auto dN_dX       = this->shape_derivatives_reference(reference_coords, r, s, t, det0);
        // deformation gradient at current position (3 x 3)
        const auto F           = this->deformation_gradient(reference_coords, current_coords, r, s, t);
        // strain displacement relationship (6 x 3N)
        const auto B           = this->green_lagrange_strain_displacement(dN_dX, F);
        // getting the current 2nd piola kirchhoff stress
        const auto S_voigt     = ip_stress.row_vec6(row);
        // getting local force of a volume element (B^T S)
        const auto local_force = B.transpose() * S_voigt * (det0 * w);

        for (Index a = 0; a < N; ++a) {
            const Index node_id = static_cast<Index>(node_ids[a]);
            for (Index d = 0; d < D; ++d) {
                node_forces(node_id, d) += local_force(D * a + d);
            }
        }
    }
}

template<Index N>
bool SolidElement<N>::compute_shear_flow(Field& shear_flow,
                                         const Field& displacement,
                                         int offset) {
    (void) shear_flow;
    (void) displacement;
    (void) offset;
    return false;
}

template<Index N>
bool SolidElement<N>::compute_beam_section_forces(Field& section_forces,
                                                  const Field& displacement,
                                                  int offset) {
    (void) section_forces;
    (void) displacement;
    (void) offset;
    return false;
}

template<Index N>
bool SolidElement<N>::compute_shell_section_forces(Field& section_forces,
                                                   Field& contribution_count,
                                                   const Field& displacement) {
    (void) section_forces;
    (void) contribution_count;
    (void) displacement;
    return false;
}

//-----------------------------------------------------------------------------
// compute_compliance
//-----------------------------------------------------------------------------
template<Index N>
void
SolidElement<N>::compute_compliance(Field& displacement, Field& result) {
    Precision buffer[D * N * D * N];
    auto K = stiffness(buffer);

    auto local_disp_mat = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    Precision strain_energy = local_displacement.dot((K * local_displacement));
    result(elem_id, 0) = strain_energy;
}

/**
 * Goal is to compute the derivative of the compliance from a linear solution w.r.t the three angles which classify
 * the additional rotation going into the section. The derivation with ()' = del () / del (alpha_i)
 * Using:
 *              K u = f
 *      K' u + K u' = 0 (chain rule)
 *
 * And the definition of compliance:
 *              J  = f^T u
 *              J' = f^T u'                  (inserting u' from above)
 *                 = f^T (-K^-1 K' u)        (inserting f^T = (Ku)^T = u^T K^T = u^T K (since K = K^T)
 *                 = - u^T K K^-1 K' u       (K K^-1 = I)
 *                 = - u^T K' u
 *
 * Using the definition of stiffness:
 *              K  = ∫B^T C_tan B dV
 *
 *              J' = - u^T ∫B^T C'_tan B dV u
 *                 = - ∫ u^T B^T C'_tan B u dV
 *                 = - ∫ eps^T C'_tan eps dV
*/
template<Index N>
void SolidElement<N>::compute_compliance_angle_derivative(Field& displacement, Field& result) {
    if (!this->_model_data || !this->_model_data->material_orientation) {
        return;
    }

    auto angles_field = this->_model_data->material_orientation;
    logging::error(angles_field->components == 3,
        "Field '", angles_field->name, "': material orientation requires 3 components");

    const Index row    = static_cast<Index>(this->elem_id);
    const Vec3  angles = angles_field->row_vec3(row);

    const Mat3 additional_rotation = cos::RectangularSystem::euler(
        angles(0),
        angles(1),
        angles(2)
    ).get_axes(Vec3::Zero());

    const std::array<Mat3, 3> additional_rotation_derivatives {
        cos::RectangularSystem::derivative_rot_x(angles(0), angles(1), angles(2)),
        cos::RectangularSystem::derivative_rot_y(angles(0), angles(1), angles(2)),
        cos::RectangularSystem::derivative_rot_z(angles(0), angles(1), angles(2))
    };

    // Initialize element state
    auto local_disp_mat     = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    const auto reference_coords = this->node_coords_reference();
    const auto current_coords   = this->node_coords_current();
    const auto& scheme          = this->integration_scheme();

    const Precision scaling    = element_stiffness_scale();
    Vec3            derivative = Vec3::Zero();

    // Integrate compliance derivatives
    for (Index n = 0; n < scheme.count(); ++n) {
        const Precision r = scheme.get_point(n).r;
        const Precision s = scheme.get_point(n).s;
        const Precision t = scheme.get_point(n).t;
        const Precision w = scheme.get_point(n).w;
        Precision det;

        // compute strain displacement matrix B
        const StaticMatrix<n_strain, D * N> B = this->strain_displacements( current_coords, r, s, t, det);
        // compute (small) strains
        const StaticVector<n_strain> strain   = B * local_displacement;

        // get the position of the point in reference coordinates, needed for the transformation
        // inside the sections coordinate system
        const Vec3 position_reference = this->interpolate<D>(reference_coords, r, s, t);

        // get the tangent rotation derivatives dC/d_alpha
        const auto tangent_derivatives = get_section()->tangent_rotation_derivatives(
            position_reference,
            additional_rotation,
            additional_rotation_derivatives
        );

        // store in the final derivative
        for (Index i = 0; i < 3; ++i) {
            derivative(i) += scaling * w * strain.dot(tangent_derivatives[i] * strain) * det;
        }
    }

    result(elem_id, 0) = derivative(0);
    result(elem_id, 1) = derivative(1);
    result(elem_id, 2) = derivative(2);
}

}  // namespace fem::model
