/**
 * @file beam.inl
 * @brief Implements the common templated base for structural beam elements.
 *
 * @see beam.h
 */

namespace fem {
namespace model {

template<Index N>
BeamElement<N>::BeamElement(ID elem_id, std::array<ID, N> node_ids_in)
    : StructuralElement(elem_id), node_ids(node_ids_in) {}

template<Index N>
BeamElement<N>::~BeamElement() = default;

template<Index N>
BeamSection* BeamElement<N>::get_section() {
    logging::error(this->_section != nullptr,
        "Section not set for element ", this->elem_id);
    logging::error(this->_section->template as<BeamSection>() != nullptr,
        "Section is not a beam section for element ", this->elem_id);
    return this->_section->template as<BeamSection>();
}

template<Index N>
Profile* BeamElement<N>::get_profile() {
    return get_section()->profile_.get();
}

template<Index N>
material::MaterialPtr BeamElement<N>::get_material() {
    BeamSection* section = get_section();
    logging::error(section->material_ != nullptr,
        "Material not set for element ", this->elem_id);
    return section->material_;
}

template<Index N>
material::IsotropicElasticity* BeamElement<N>::get_elasticity() {
    BeamSection* section = get_section();
    logging::error(section->material_ != nullptr,
        "Material not set for element ", this->elem_id);
    logging::error(section->material_->has_elasticity(),
        "Material has no elasticity assigned");
    logging::error(section->material_->elasticity()->template as<material::IsotropicElasticity>() != nullptr,
        "Material is not isotropic for element ", this->elem_id);
    return section->material_->elasticity()->template as<material::IsotropicElasticity>();
}

template<Index N>
Vec3 BeamElement<N>::orientation_direction() {
    constexpr Precision kOrientationEps = static_cast<Precision>(1e-12);
    BeamSection* section = get_section();
    logging::error(section->direction_.norm() > kOrientationEps,
        "Beam element ", this->elem_id, " requires a non-zero section n1 vector");
    return section->direction_.normalized();
}

template<Index N>
Precision BeamElement<N>::length() {
    Precision l = 0;
    for (Index i = 0; i < N - 1; ++i) {
        const Vec3 x0 = this->node_position(static_cast<ID>(i));
        const Vec3 x1 = this->node_position(static_cast<ID>(i + 1));
        l += (x0 - x1).norm();
    }
    return l;
}

template<Index N>
Precision BeamElement<N>::volume() {
    return get_profile()->area_ * length();
}

template<Index N>
Mat3 BeamElement<N>::base_rotation_matrix() {
    Vec3 x = (this->node_position(1) - this->node_position(0)).normalized();
    Vec3 y = orientation_direction();
    Vec3 z = x.cross(y).normalized();
    y = z.cross(x).normalized();

    Mat3 mat{};
    mat(0, 0) = x(0); mat(0, 1) = x(1); mat(0, 2) = x(2);
    mat(1, 0) = y(0); mat(1, 1) = y(1); mat(1, 2) = y(2);
    mat(2, 0) = z(0); mat(2, 1) = z(1); mat(2, 2) = z(2);
    return mat;
}

template<Index N>
Precision BeamElement<N>::principal_angle() {
    Profile* pr = get_profile();
    const Precision Iy  = pr->inertia_y_;
    const Precision Iz  = pr->inertia_z_;
    const Precision Iyz = pr->product_inertia_yz_;

    const Precision scale = std::max<Precision>(Precision(1), std::abs(Iy) + std::abs(Iz));
    if (std::abs(Iyz) <= scale * Precision(1e-14)) return Precision(0);
    return Precision(0.5) * std::atan2(Precision(2) * Iyz, Iz - Iy);
}

template<Index N>
Mat3 BeamElement<N>::principal_rotation_matrix() {
    Mat3 Rb = base_rotation_matrix();
    const Precision phi = principal_angle();
    if (phi == Precision(0)) return Rb;

    Eigen::Matrix<Precision, 1, 3> rx = Rb.row(0);
    Eigen::Matrix<Precision, 1, 3> ry = Rb.row(1);
    Eigen::Matrix<Precision, 1, 3> rz = Rb.row(2);

    const Precision c = std::cos(phi);
    const Precision s = std::sin(phi);

    Eigen::Matrix<Precision, 1, 3> ry_p =  c * ry + s * rz;
    Eigen::Matrix<Precision, 1, 3> rz_p = -s * ry + c * rz;

    Mat3 Rp = Rb;
    Rp.row(0) = rx;
    Rp.row(1) = ry_p;
    Rp.row(2) = rz_p;
    return Rp;
}

template<Index N>
Mat3 BeamElement<N>::rotation_matrix() {
    return principal_rotation_matrix();
}

template<Index N>
StaticMatrix<6, 6> BeamElement<N>::rigid_offset_6(Precision dy, Precision dz) {
    StaticMatrix<6, 6> B = StaticMatrix<6, 6>::Identity();
    B(0, 4) =  dz;
    B(0, 5) = -dy;
    B(1, 3) = -dz;
    B(2, 3) =  dy;
    return B;
}

template<Index N>
StaticMatrix<N * 6, N * 6> BeamElement<N>::rigid_offset_N(Precision dy, Precision dz) {
    StaticMatrix<N * 6, N * 6> B = StaticMatrix<N * 6, N * 6>::Identity();
    const StaticMatrix<6, 6> Bi = rigid_offset_6(dy, dz);
    for (Index a = 0; a < N; ++a) {
        B.template block<6, 6>(a * 6, a * 6) = Bi;
    }
    return B;
}

template<Index N>
void BeamElement<N>::rotate_yz_to_principal(Precision phi, Precision& y, Precision& z) {
    if (phi == Precision(0)) return;
    const Precision c = std::cos(phi);
    const Precision s = std::sin(phi);
    const Precision y_p =  c * y + s * z;
    const Precision z_p = -s * y + c * z;
    y = y_p;
    z = z_p;
}

template<Index N>
StaticMatrix<N * 6, N * 6> BeamElement<N>::transformation() {
    StaticMatrix<N * 6, N * 6> T;
    T.setZero();
    Mat3 R = rotation_matrix();

    for (Index i = 0; i < N; ++i) {
        for (Dim j = 0; j < 3; ++j) {
            for (Dim k = 0; k < 3; ++k) {
                T(i * 6 + j,     i * 6 + k)     = R(j, k);
                T(i * 6 + j + 3, i * 6 + k + 3) = R(j, k);
            }
        }
    }
    return T;
}

template<Index N>
StaticMatrix<N * 6, N * 6> BeamElement<N>::transformation_base() {
    StaticMatrix<N * 6, N * 6> T;
    T.setZero();
    Mat3 R = base_rotation_matrix();

    for (Index i = 0; i < N; ++i) {
        for (Dim j = 0; j < 3; ++j) {
            for (Dim k = 0; k < 3; ++k) {
                T(i * 6 + j,     i * 6 + k)     = R(j, k);
                T(i * 6 + j + 3, i * 6 + k + 3) = R(j, k);
            }
        }
    }
    return T;
}

template<Index N>
MapMatrix BeamElement<N>::stiffness(Precision* buffer) {
    MapMatrix result(buffer, N * 6, N * 6);
    result = stiffness_impl();
    return result;
}

template<Index N>
MapMatrix BeamElement<N>::stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) {
    MapMatrix result(buffer, N * 6, N * 6);
    result = stiffness_geom_impl(ip_stress, ip_start_idx);
    return result;
}

template<Index N>
MapMatrix BeamElement<N>::mass(Precision* buffer) {
    MapMatrix result(buffer, N * 6, N * 6);
    result = mass_impl();
    return result;
}

template<Index N>
void BeamElement<N>::compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) {
    (void) node_forces;
    (void) ip_stress;
    logging::error(false,
        "BeamElement: compute_internal_force_nonlinear is not implemented yet for element ", this->elem_id);
}

template<Index N>
void BeamElement<N>::compute_stress_state(
    Field&       stress_state,
    const Field& displacement,
    int          offset,
    bool         use_green_lagrange_nl
) {
    RowMatrix rst = stress_strain_ip_rst();
    if (rst.rows() == 0) return;
    compute_stress_strain(nullptr, &stress_state, displacement, rst, offset, use_green_lagrange_nl);
}

template<Index N>
void BeamElement<N>::apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) {
    (void) node_loads;
    (void) node_temp;
    (void) ref_temp;
}

template<Index N>
ElDofs BeamElement<N>::dofs() const {
    return ElDofs{true, true, true, true, true, true};
}

template<Index N>
Dim BeamElement<N>::dimensions() const {
    return 3;
}

template<Index N>
Dim BeamElement<N>::n_nodes() const {
    return N;
}

template<Index N>
Dim BeamElement<N>::num_ip() const {
    return 1;
}

template<Index N>
const ID* BeamElement<N>::nodes() const {
    return node_ids.data();
}

template<Index N>
SurfacePtr BeamElement<N>::surface(ID surface_id) {
    (void) surface_id;
    return nullptr;
}

template<Index N>
bool BeamElement<N>::compute_beam_section_forces(
    Field&       section_forces,
    const Field& displacement,
    int          offset
) {
    Eigen::Matrix<Precision, N * 6, 1> u_global;
    for (Index i = 0; i < N; ++i) {
        const ID nid = node_ids[i];
        const Vec6 row = displacement.row_vec6(static_cast<Index>(nid));
        for (Index d = 0; d < 6; ++d) {
            u_global(i * 6 + d) = row(d);
        }
    }

    const auto K_global = stiffness_impl();
    const auto T_out    = transformation_base();
    const auto f_global = K_global * u_global;
    const auto q_local  = T_out * f_global;

    for (Index i = 0; i < N; ++i) {
        for (Index d = 0; d < 6; ++d) {
            section_forces(static_cast<Index>(offset) + i, d)
                = q_local(i * 6 + d) * ((i == 1 && N == 2) ? 1 : -1);
        }
    }
    return true;
}

template<Index N>
Precision BeamElement<N>::integrate_scalar_field(bool scale_by_density, const ScalarField& field) {
    const Precision L = length();
    const Precision A = get_profile()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return Precision(0);

    Vec3 x_mid = Vec3::Zero();
    for (Index i = 0; i < N; ++i) x_mid += this->node_position(static_cast<ID>(i));
    x_mid /= static_cast<Precision>(N);

    Precision rho = Precision(1);
    if (scale_by_density) {
        auto mat = get_material();
        logging::error(mat && mat->has_density(),
            "BeamElement: material density is required when scale_by_density=true for element ", this->elem_id);
        rho = mat->get_density();
    }
    return field(x_mid) * (rho * A * L);
}

template<Index N>
Vec3 BeamElement<N>::integrate_vector_field(bool scale_by_density, const VecField& field) {
    const Precision L = length();
    const Precision A = get_profile()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return Vec3::Zero();

    Vec3 x_mid = Vec3::Zero();
    for (Index i = 0; i < N; ++i) x_mid += this->node_position(static_cast<ID>(i));
    x_mid /= static_cast<Precision>(N);

    Precision rho = Precision(1);
    if (scale_by_density) {
        auto mat = get_material();
        logging::error(mat && mat->has_density(),
            "BeamElement: material density is required when scale_by_density=true for element ", this->elem_id);
        rho = mat->get_density();
    }
    return field(x_mid) * (rho * A * L);
}

template<Index N>
void BeamElement<N>::integrate_vector_field(
    Field&          node_loads,
    bool            scale_by_density,
    const VecField& field
) {
    const Precision L = length();
    const Precision A = get_profile()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return;

    Vec3 x_mid = Vec3::Zero();
    for (Index i = 0; i < N; ++i) x_mid += this->node_position(static_cast<ID>(i));
    x_mid /= static_cast<Precision>(N);

    Precision rho = Precision(1);
    if (scale_by_density) {
        auto mat = get_material();
        logging::error(mat && mat->has_density(),
            "BeamElement: material density is required when scale_by_density=true for element ", this->elem_id);
        rho = mat->get_density();
    }

    const Vec3 F = field(x_mid) * (rho * A * L);
    const Precision share = Precision(1) / static_cast<Precision>(N);
    for (Index i = 0; i < N; ++i) {
        const ID n_id = node_ids[i];
        node_loads(n_id, 0) += share * F(0);
        node_loads(n_id, 1) += share * F(1);
        node_loads(n_id, 2) += share * F(2);
    }
}

template<Index N>
Mat3 BeamElement<N>::integrate_tensor_field(bool scale_by_density, const TenField& field) {
    const Precision L = length();
    const Precision A = get_profile()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return Mat3::Zero();

    Vec3 x_mid = Vec3::Zero();
    for (Index i = 0; i < N; ++i) x_mid += this->node_position(static_cast<ID>(i));
    x_mid /= static_cast<Precision>(N);

    Precision rho = Precision(1);
    if (scale_by_density) {
        auto mat = get_material();
        logging::error(mat && mat->has_density(),
            "BeamElement: material density is required when scale_by_density=true for element ", this->elem_id);
        rho = mat->get_density();
    }
    return field(x_mid) * (rho * A * L);
}

} // namespace model
} // namespace fem
