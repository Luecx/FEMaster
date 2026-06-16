/**
 * @file element_solid.ipp
 * @brief Implements the SolidElement class template.
 *
 * This file contains element-level routines for solid elements, including
 * stress/strain recovery, internal force computation, compliance evaluation,
 * and topology-orientation derivatives.
 *
 * @date Created on 12.06.2023
 * @author Finn Eggers
 */

#pragma once

namespace fem::model {

namespace detail_solid_nonlinear {

/**
 * @brief Converts a stress-like Voigt vector to a symmetric 3x3 tensor.
 *
 * Voigt order:
 *
 *     [11, 22, 33, 23, 31, 12]
 */
inline Mat3 voigt_to_tensor(const Vec6& v) {
    Mat3 tensor;

    tensor << v(0), v(5), v(4),
              v(5), v(1), v(3),
              v(4), v(3), v(2);

    return tensor;
}

/**
 * @brief Converts a symmetric 3x3 tensor to a stress-like Voigt vector.
 *
 * Voigt order:
 *
 *     [11, 22, 33, 23, 31, 12]
 */
inline Vec6 tensor_to_voigt(const Mat3& tensor) {
    Vec6 v;

    v << tensor(0, 0),
         tensor(1, 1),
         tensor(2, 2),
         tensor(1, 2),
         tensor(2, 0),
         tensor(0, 1);

    return v;
}

/**
 * @brief Converts the Green-Lagrange strain tensor to engineering Voigt strain.
 *
 * Voigt order:
 *
 *     [E11, E22, E33, 2E23, 2E31, 2E12]
 */
inline Vec6 green_lagrange_to_voigt(const Mat3& strain) {
    Vec6 v;

    v << strain(0, 0),
         strain(1, 1),
         strain(2, 2),
         Precision(2) * strain(1, 2),
         Precision(2) * strain(2, 0),
         Precision(2) * strain(0, 1);

    return v;
}

/**
 * @brief Computes the displacement gradient with respect to reference coordinates.
 */
template<Index N>
inline Mat3 displacement_gradient(const StaticMatrix<N, 3>& displacement,
                                  const StaticMatrix<N, 3>& dN_dX) {
    Mat3 grad_u = Mat3::Zero();

    for (Index a = 0; a < N; ++a) {
        grad_u += displacement.row(a).transpose() * dN_dX.row(a);
    }

    return grad_u;
}

/**
 * @brief Computes the Cauchy stress from Green-Lagrange strain using St. Venant-Kirchhoff kinematics.
 */
inline Vec6 compute_cauchy_stress(const Vec6& material_strain,
                                  const Mat3& F,
                                  const Mat6& C) {
    const Precision J             = F.determinant();
    const Vec6      second_pk_vec = C * material_strain;
    const Mat3      second_pk     = voigt_to_tensor(second_pk_vec);
    const Mat3      cauchy        = (F * second_pk * F.transpose()) / J;

    return tensor_to_voigt(cauchy);
}

} // namespace detail_solid_nonlinear


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

    const auto node_coords        = this->node_coords_global();
    const auto local_displacement = this->nodal_data<3>(displacement);

    const auto local_disp_mat     = StaticMatrix<3, N>(local_displacement.transpose());
    const auto local_disp_vec     = Eigen::Map<const StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    for (Index n = 0; n < rst.rows(); ++n) {
        const Precision r   = rst(n, 0);
        const Precision s   = rst(n, 1);
        const Precision t   = rst(n, 2);
        const Index     row = static_cast<Index>(offset) + n;

        if (!use_green_lagrange_nl) {
            Precision det = Precision(0);

            const StaticMatrix<n_strain, D * N> B =
                this->strain_displacements(node_coords, r, s, t, det, false);

            if (det <= Precision(0) || std::isnan(det) || std::isinf(det)) {
                continue;
            }

            const Vec6 global_strain = B * local_disp_vec;

            Vec6 local_strain;
            Vec6 local_stress;

            material_stress_strain(r, s, t, global_strain, local_stress, local_strain);

            for (Dim j = 0; j < n_strain; ++j) {
                if (strain != nullptr) {
                    (*strain)(row, j) = local_strain(j);
                }

                if (stress != nullptr) {
                    (*stress)(row, j) = local_stress(j);
                }
            }

            continue;
        }

        const StaticMatrix<N, D> dN_local = this->shape_derivative(r, s, t);
        const StaticMatrix<D, D> J0       = this->jacobian(node_coords, r, s, t);
        const Precision          detJ0    = J0.determinant();

        logging::error(detJ0 > Precision(0),
                       "negative reference determinant encountered in nonlinear stress for element ", elem_id,
                       "\ndet        : ", detJ0,
                       "\nCoordinates: ", node_coords);

        const StaticMatrix<N, D> dN_dX =
            (J0.inverse() * dN_local.transpose()).transpose();

        const Mat3 grad_u = detail_solid_nonlinear::displacement_gradient<N>(
            local_displacement,
            dN_dX
        );

        const Mat3      F = Mat3::Identity() + grad_u;
        const Precision J = F.determinant();

        logging::error(J > Precision(0),
                       "non-positive deformation gradient determinant in element ", elem_id,
                       "\nJ: ", J,
                       "\nF: ", F);

        const Mat3 green_lagrange =
            Precision(0.5) * (F.transpose() * F - Mat3::Identity());

        const Vec6 green_voigt =
            detail_solid_nonlinear::green_lagrange_to_voigt(green_lagrange);

        const Vec6 cauchy_voigt =
            detail_solid_nonlinear::compute_cauchy_stress(
                green_voigt,
                F,
                material_matrix(r, s, t)
            );

        for (Dim j = 0; j < n_strain; ++j) {
            if (strain != nullptr) {
                (*strain)(row, j) = green_voigt(j);
            }

            if (stress != nullptr) {
                (*stress)(row, j) = cauchy_voigt(j);
            }
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

    compute_stress_strain(
        nullptr,
        &stress_state,
        displacement,
        rst,
        offset,
        use_green_lagrange_nl
    );
}


template<Index N>
void SolidElement<N>::compute_internal_force_nonlinear(Field& node_forces,
                                                       const Field& displacement,
                                                       const Field& ip_stress,
                                                       int ip_offset) {
    (void) displacement;

    const auto current_coords = this->node_coords_global();
    const auto scheme         = this->integration_scheme();

    for (Index n = 0; n < scheme.count(); ++n) {
        const auto&     point = scheme.get_point(n);
        const Precision r     = point.r;
        const Precision s     = point.s;
        const Precision t     = point.t;
        const Precision w     = point.w;
        const Index     row   = static_cast<Index>(ip_offset) + n;

        Precision det = Precision(0);

        const StaticMatrix<n_strain, D * N> B =
            this->strain_displacements(current_coords, r, s, t, det);

        const Vec6                sigma       = ip_stress.row_vec6(row);
        const StaticVector<D * N> local_force = B.transpose() * sigma * (det * w);

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


template<Index N>
void SolidElement<N>::compute_compliance(Field& displacement,
                                         Field& result) {
    Precision buffer[D * N * D * N];

    const auto K = stiffness(buffer);

    const auto local_disp_mat =
        StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());

    const auto local_displacement =
        Eigen::Map<const StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    const Precision strain_energy =
        local_displacement.dot(K * local_displacement);

    result(elem_id, 0) = strain_energy;
}


template<Index N>
void SolidElement<N>::compute_compliance_angle_derivative(Field& displacement,
                                                          Field& result) {
    if (!this->has_topo_orientation()) {
        return;
    }

    const auto local_disp_mat =
        StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());

    const auto local_displacement =
        Eigen::Map<const StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    const auto node_coords = this->node_coords_global();
    const auto scheme      = this->integration_scheme();

    Vec3 derivative = Vec3::Zero();

    for (Index n = 0; n < scheme.count(); ++n) {
        const auto&     point = scheme.get_point(n);
        const Precision r     = point.r;
        const Precision s     = point.s;
        const Precision t     = point.t;
        const Precision w     = point.w;

        Precision det = Precision(0);

        const StaticMatrix<n_strain, D * N> B =
            this->strain_displacements(node_coords, r, s, t, det);

        const StaticVector<n_strain> strain = B * local_displacement;
        const Mat3                   R      = this->material_basis(r, s, t);

        for (Index angle = 0; angle < 3; ++angle) {
            const Mat3 dR = this->material_basis_derivative(r, s, t, angle);

            const auto dC =
                this->material()
                    ->elasticity()
                    ->template get_transformed_derivative<3>(R, dR);

            derivative(angle) += w * strain.dot(dC * strain) * det;
        }
    }

    result(elem_id, 0) = derivative(0);
    result(elem_id, 1) = derivative(1);
    result(elem_id, 2) = derivative(2);
}

} // namespace fem::model
