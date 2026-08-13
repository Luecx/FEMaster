/**
 * @file c3d8r_base.ipp
 * @brief Basic C3D8R interpolation and reference geometry helpers.
 *
 * @author Finn Eggers
 * @date 13.08.2026
 */

C3D8R::C3D8R(ID elem_id, const std::array<ID, N>& node_ids)
    : C3D8(elem_id, node_ids) {}

std::string C3D8R::type_name() const {
    return "C3D8R";
}

const math::quadrature::Quadrature& C3D8R::integration_scheme_stiffness() const {
    static const math::quadrature::Quadrature quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_CONSTANT
    };
    return quadrature;
}

RowMatrix C3D8R::stress_strain_nodal_rst() {
    return RowMatrix::Zero(N, D);
}

C3D8R::HourglassModes C3D8R::primitive_hourglass_modes() {
    const auto local = node_coords_local();
    HourglassModes modes = HourglassModes::Zero();

    for (Index a = 0; a < N; ++a) {
        const Precision r = local(a, 0);
        const Precision s = local(a, 1);
        const Precision t = local(a, 2);
        modes(a, 0) = s * t;
        modes(a, 1) = t * r;
        modes(a, 2) = r * s;
        modes(a, 3) = r * s * t;
    }
    return modes;
}

C3D8R::GradientMatrix C3D8R::mean_reference_gradient() {
    const auto reference = node_coords_reference();
    static const math::quadrature::Quadrature full{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_QUADRATIC
    };

    GradientMatrix integrated = GradientMatrix::Zero();
    Precision volume = Precision(0);

    for (Index q = 0; q < full.count(); ++q) {
        const auto point = full.get_point(q);
        Precision det0 = Precision(0);
        const auto gradient = shape_derivatives_reference(
            reference, point.r, point.s, point.t, det0);

        logging::error(std::isfinite(det0) && det0 > Precision(0),
            "C3D8R: invalid reference determinant in element ", elem_id,
            "\ndet(J0): ", det0);

        const Precision measure = det0 * point.w;
        integrated += gradient * measure;
        volume += measure;
    }

    logging::error(std::isfinite(volume) && volume > Precision(0),
        "C3D8R: invalid reference volume in element ", elem_id,
        "\nvolume: ", volume);
    return integrated / volume;
}

Mat3 C3D8R::hourglass_reference_frame() {
    const auto reference = node_coords_reference();
    const Mat3 J0 = jacobian(reference, Precision(0), Precision(0), Precision(0));

    logging::error(std::isfinite(J0.determinant()) && J0.determinant() > Precision(0),
        "C3D8R: invalid center reference Jacobian in element ", elem_id);

    Vec3 a1 = J0.row(0).transpose();
    Vec3 a2 = J0.row(1).transpose();
    const Vec3 a3 = J0.row(2).transpose();

    logging::error(a1.norm() > Precision(0),
        "C3D8R: degenerate first reference axis in element ", elem_id);
    const Vec3 e1 = a1.normalized();

    a2 -= e1 * e1.dot(a2);
    logging::error(a2.norm() > Precision(0),
        "C3D8R: degenerate second reference axis in element ", elem_id);
    Vec3 e2 = a2.normalized();
    Vec3 e3 = e1.cross(e2).normalized();

    if (e3.dot(a3) < Precision(0)) {
        e2 = -e2;
        e3 = -e3;
    }

    Mat3 frame;
    frame.col(0) = e1;
    frame.col(1) = e2;
    frame.col(2) = e3;

    logging::error(frame.allFinite() && frame.determinant() > Precision(0),
        "C3D8R: invalid hourglass reference frame in element ", elem_id);
    return frame;
}

Mat3 C3D8R::hourglass_geometry_integrals(const Mat3& frame) {
    const auto reference = node_coords_reference();
    const auto local = node_coords_local();
    const GradientMatrix reference_local = reference * frame;

    StaticVector<D> L = StaticVector<D>::Zero();
    for (Dim i = 0; i < D; ++i) {
        L(i) = local.col(i).dot(reference_local.col(i));
        logging::error(std::isfinite(L(i)) && L(i) > Precision(0),
            "C3D8R: invalid physical hourglass length in element ", elem_id);
    }

    Mat3 H = Mat3::Zero();
    for (Dim i = 0; i < D; ++i) {
        const Dim j = (i + 1) % D;
        const Dim k = (i + 2) % D;
        H(i, i) = L(j) * L(k) / (Precision(3) * L(i));
    }

    for (Dim i = 0; i < D; ++i) {
        for (Dim j = i + 1; j < D; ++j) {
            const Dim k = D - i - j;
            H(i, j) = H(j, i) = L(k) / Precision(3);
        }
    }
    return H;
}
