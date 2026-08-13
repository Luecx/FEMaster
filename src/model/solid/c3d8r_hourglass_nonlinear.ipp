/**
 * @file c3d8r_hourglass_nonlinear.ipp
 * @brief Evaluates the objective Total-Lagrangian hourglass potential.
 *
 * @author Finn Eggers
 * @date 13.08.2026
 */

/**
 * Evaluates
 *
 *     Pi_hg = 1/2 z^T A_hg z,
 *     z_alpha = R0^T F_bar^T q_alpha,
 *
 * with `q_alpha = sum gamma_a_alpha u_a`. The pull-back is invariant under a
 * superposed rigid rotation. Residual and tangent are exact first and second
 * derivatives of the same potential.
 */
void C3D8R::hourglass_response(Vector24& force, Matrix24* tangent) {
    build_hourglass_reference_data();

    const GradientMatrix displacement =
        node_coords_current() - node_coords_reference();
    const Mat3 F =
        Mat3::Identity() + displacement.transpose() * mean_gradient_cache;

    logging::error(F.allFinite() && F.determinant() > Precision(0),
        "C3D8R: non-positive mean deformation gradient in element ", elem_id,
        "\ndet(F_bar): ", F.determinant());

    const StaticMatrix<D, n_hg_mode> q =
        displacement.transpose() * hourglass_modes_cache;
    const StaticMatrix<D, n_hg_mode> z_matrix =
        hourglass_frame_cache.transpose() * F.transpose() * q;

    ModalVector z = ModalVector::Zero();
    for (Dim i = 0; i < D; ++i) {
        for (Index alpha = 0; alpha < n_hg_mode; ++alpha) {
            z(n_hg_mode * i + alpha) = z_matrix(i, alpha);
        }
    }

    const ModalVector p = modal_stiffness_cache * z;
    ModalJacobian J = ModalJacobian::Zero();

    // Exact first derivative dz/du.
    for (Index a = 0; a < N; ++a) {
        const Vec3 g = mean_gradient_cache.row(a).transpose();

        for (Index alpha = 0; alpha < n_hg_mode; ++alpha) {
            const Mat3 dz_du = hourglass_frame_cache.transpose() * (
                g * q.col(alpha).transpose() +
                hourglass_modes_cache(a, alpha) * F.transpose()
            );

            for (Dim i = 0; i < D; ++i) {
                for (Dim d = 0; d < D; ++d) {
                    J(n_hg_mode * i + alpha, D * a + d) = dz_du(i, d);
                }
            }
        }
    }

    force = J.transpose() * p;
    logging::error(force.allFinite(),
        "C3D8R: non-finite finite-strain hourglass force in element ", elem_id);

    if (tangent == nullptr) return;

    *tangent = J.transpose() * modal_stiffness_cache * J;

    // Contract the exact second derivative of z = R0^T F_bar^T q with p.
    for (Index a = 0; a < N; ++a) {
        const Vec3 ga = mean_gradient_cache.row(a).transpose();

        for (Index b = 0; b < N; ++b) {
            const Vec3 gb = mean_gradient_cache.row(b).transpose();
            Precision scalar = Precision(0);

            for (Index alpha = 0; alpha < n_hg_mode; ++alpha) {
                const Vec3 h = hourglass_frame_cache.transpose() * (
                    ga * hourglass_modes_cache(b, alpha) +
                    gb * hourglass_modes_cache(a, alpha)
                );

                for (Dim i = 0; i < D; ++i) {
                    scalar += p(n_hg_mode * i + alpha) * h(i);
                }
            }

            for (Dim d = 0; d < D; ++d) {
                (*tangent)(D * a + d, D * b + d) += scalar;
            }
        }
    }

    *tangent = Precision(0.5) * (*tangent + tangent->transpose());
    logging::error(tangent->allFinite(),
        "C3D8R: non-finite finite-strain hourglass tangent in element ", elem_id);
}
