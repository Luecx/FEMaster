/**
 * @file c3d8r_hourglass_material.ipp
 * @brief Extracts state-neutral physical hourglass material parameters.
 *
 * @author Finn Eggers
 * @date 13.08.2026
 */

StaticVector<2> C3D8R::hourglass_material_parameters() {
    const Index row = mp_index(0);
    Field& field = *_model_data->material_state;

    std::vector<Precision> saved(static_cast<std::size_t>(field.components));
    for (Index c = 0; c < field.components; ++c) {
        saved[static_cast<std::size_t>(c)] = field(row, c);
    }

    Precision* state = &field(row, 0);
    const Mat6 C = material_tangent_reference(
        Precision(0), Precision(0), Precision(0), state);

    for (Index c = 0; c < field.components; ++c) {
        field(row, c) = saved[static_cast<std::size_t>(c)];
    }

    const Precision mu =
        (C(3, 3) + C(4, 4) + C(5, 5)) / Precision(3);
    const Precision lambda =
        (C(0, 1) + C(1, 0) + C(0, 2) + C(2, 0) + C(1, 2) + C(2, 1)) /
        Precision(6);
    const Precision nu = lambda / (Precision(2) * (lambda + mu));

    logging::error(std::isfinite(mu) && mu > Precision(0),
        "C3D8R: invalid initial shear stiffness in element ", elem_id);
    logging::error(std::isfinite(nu) && nu > Precision(-1) && nu < Precision(0.5),
        "C3D8R: invalid isotropic-equivalent Poisson ratio in element ", elem_id);

    StaticVector<2> result;
    result << mu, nu;
    return result;
}
