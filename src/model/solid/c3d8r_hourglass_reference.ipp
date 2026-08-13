/**
 * @file c3d8r_hourglass_reference.ipp
 * @brief Builds the reference Belytschko-Bindeman modal stiffness.
 *
 * @author Finn Eggers
 * @date 13.08.2026
 */

void C3D8R::build_hourglass_reference_data() {
    if (hourglass_reference_cached) return;

    const auto reference = node_coords_reference();
    mean_gradient_cache = mean_reference_gradient();
    hourglass_modes_cache =
        (StaticMatrix<N, N>::Identity() - mean_gradient_cache * reference.transpose()) *
        primitive_hourglass_modes();
    hourglass_frame_cache = hourglass_reference_frame();

    const Mat3 H = hourglass_geometry_integrals(hourglass_frame_cache);
    const auto material = hourglass_material_parameters();
    const Precision mu = material(0);
    const Precision nu = material(1);
    const Precision nf = (Precision(1) + nu) / (Precision(1) - nu);
    const Precision pf = nu / (Precision(1) - nu);

    modal_stiffness_cache.setZero();
    for (Dim i = 0; i < D; ++i) {
        const Dim j = (i + 1) % D;
        const Dim k = (i + 2) % D;
        modal_stiffness_cache(n_hg_mode * i + j, n_hg_mode * i + j) += H(i, i) * nf;
        modal_stiffness_cache(n_hg_mode * i + k, n_hg_mode * i + k) += H(i, i) * nf;
        modal_stiffness_cache(n_hg_mode * i + 3, n_hg_mode * i + 3) += H(i, i) / Precision(3);
        modal_stiffness_cache(n_hg_mode * i + i, n_hg_mode * i + i) +=
            Precision(0.5) * (H(j, j) + H(k, k));
    }

    for (Dim i = 0; i < D; ++i) {
        for (Dim j = 0; j < D; ++j) {
            if (i == j) continue;
            modal_stiffness_cache(n_hg_mode * i + j, n_hg_mode * j + i) += H(i, j) * pf;
            modal_stiffness_cache(n_hg_mode * i + i, n_hg_mode * j + j) +=
                Precision(0.5) * H(i, j);
        }
    }

    modal_stiffness_cache *= Precision(2) * mu;
    modal_stiffness_cache = Precision(0.5) *
        (modal_stiffness_cache + modal_stiffness_cache.transpose());

    ModalJacobian J0 = ModalJacobian::Zero();
    for (Index a = 0; a < N; ++a) {
        for (Index alpha = 0; alpha < n_hg_mode; ++alpha) {
            for (Dim i = 0; i < D; ++i) {
                for (Dim d = 0; d < D; ++d) {
                    J0(n_hg_mode * i + alpha, D * a + d) =
                        hourglass_modes_cache(a, alpha) * hourglass_frame_cache(d, i);
                }
            }
        }
    }

    hourglass_stiffness_cache = J0.transpose() * modal_stiffness_cache * J0;
    hourglass_stiffness_cache = Precision(0.5) *
        (hourglass_stiffness_cache + hourglass_stiffness_cache.transpose());
    logging::error(hourglass_stiffness_cache.allFinite(),
        "C3D8R: non-finite physical hourglass stiffness in element ", elem_id);

    hourglass_reference_cached = true;
}

C3D8R::Matrix24 C3D8R::hourglass_stiffness() {
    build_hourglass_reference_data();
    return hourglass_stiffness_cache;
}
