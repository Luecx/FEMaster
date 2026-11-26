/**
 * @file beam.inl
 * @brief Implements the common templated base for structural beam elements.
 *
 * The beam subsystem shares section frames, rigid offsets, reference-state
 * validation and linear mechanical/result operators through BeamElement<N>.
 * Explicit zero-valued base displacements and an omitted base both represent
 * the reference configuration for those linear operators. Nonzero bases remain
 * the responsibility of a finite-rotation formulation in a derived element.
 *
 * Concrete beam elements supply elastic, geometric and mass matrices. Assembly
 * and result writers own global fields, while these definitions gather and
 * scatter the connected nodal quantities without changing the base state.
 *
 * @see beam.h
 * @see BeamElement
 */

namespace fem {
namespace model {

template<Index N>
BeamElement<N>::BeamElement(Index elem_id, std::array<Index, N> node_ids_in)
    : StructuralElement(elem_id), node_ids(node_ids_in) {}

template<Index N>
BeamElement<N>::~BeamElement() = default;

template<Index N>
BeamSection* BeamElement<N>::get_section() {
    logging::error(this->_section != nullptr,
        "Section not set for element ", this->label());
    logging::error(this->_section->template as<BeamSection>() != nullptr,
        "Section is not a beam section for element ", this->label());
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
        "Material not set for element ", this->label());
    return section->material_;
}

template<Index N>
material::IsotropicElasticity* BeamElement<N>::get_elasticity() {
    BeamSection* section = get_section();
    logging::error(section->material_ != nullptr,
        "Material not set for element ", this->label());
    logging::error(section->material_->has_elasticity(),
        "Material has no elasticity assigned");
    logging::error(section->material_->elasticity()->template as<material::IsotropicElasticity>() != nullptr,
        "Material is not isotropic for element ", this->label());
    return section->material_->elasticity()->template as<material::IsotropicElasticity>();
}

template<Index N>
Vec3 BeamElement<N>::orientation_direction() {
    constexpr Precision kOrientationEps = static_cast<Precision>(1e-12);
    BeamSection* section = get_section();
    logging::error(section->direction_.norm() > kOrientationEps,
        "Beam element ", this->label(), " requires a non-zero section n1 vector");
    return section->direction_.normalized();
}

template<Index N>
Precision BeamElement<N>::length() {
    Precision l = 0;
    for (Index i = 0; i < N - 1; ++i) {
        const Vec3 x0 = this->node_position(i);
        const Vec3 x1 = this->node_position(i + 1);
        l += (x0 - x1).norm();
    }
    return l;
}

/**
 * Accumulates the centerline length over consecutive reference nodes.
 *
 * @return Reference centerline length L0.
 */
template<Index N>
Precision BeamElement<N>::length_reference() {
    logging::error(this->_model_data != nullptr && this->_model_data->positions_reference != nullptr,
        "Beam element ", this->label(), ": reference positions are not initialized");

    // Sum the straight segment lengths of the undeformed centerline
    const auto& reference = *this->_model_data->positions_reference;

    Precision l = 0;
    for (Index i = 0; i < N - 1; ++i) {
        const Vec3 x0 = reference.row_vec3(node_ids[i]);
        const Vec3 x1 = reference.row_vec3(node_ids[i + 1]);
        l += (x0 - x1).norm();
    }
    return l;
}

/**
 * Computes the arithmetic mean of the reference nodal positions.
 *
 * The one-point field integration of beams samples distributed fields at this
 * point and distributes the resultant equally to the nodes.
 *
 * @return Reference nodal midpoint X_m.
 */
template<Index N>
Vec3 BeamElement<N>::midpoint_reference() {
    logging::error(this->_model_data != nullptr && this->_model_data->positions_reference != nullptr,
        "Beam element ", this->label(), ": reference positions are not initialized");

    // Average the undeformed nodal positions
    const auto& reference = *this->_model_data->positions_reference;

    Vec3 x_mid = Vec3::Zero();
    for (Index i = 0; i < N; ++i) {
        x_mid += reference.row_vec3(node_ids[i]);
    }
    return x_mid / static_cast<Precision>(N);
}

template<Index N>
Precision BeamElement<N>::volume() {
    return get_profile()->area_ * length();
}

template<Index N>
Mat3 BeamElement<N>::base_rotation_matrix() {
    // Validate the element axis before normalizing its direction
    const Vec3 axis = this->node_position(1) - this->node_position(0);
    logging::error(axis.norm() > Precision(1e-12),
        "Beam element ", this->label(), ": degenerate beam axis");

    Vec3 x = axis.normalized();
    Vec3 y = orientation_direction();
    Vec3 z = x.cross(y);

    // Reject a section direction parallel to the element axis
    logging::error(z.norm() > Precision(1e-12),
        "Beam element ", this->label(), ": section n1 is parallel to the beam axis");

    z.normalize();
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
        for (Index j = 0; j < 3; ++j) {
            for (Index k = 0; k < 3; ++k) {
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
        for (Index j = 0; j < 3; ++j) {
            for (Index k = 0; k < 3; ++k) {
                T(i * 6 + j,     i * 6 + k)     = R(j, k);
                T(i * 6 + j + 3, i * 6 + k + 3) = R(j, k);
            }
        }
    }
    return T;
}

/**
 * Validates the undeformed base required by the shared linear beam operators.
 *
 * A null pointer denotes the implicit reference state. An explicit field is
 * equivalent only when all six global generalized displacement components at
 * every connected node are exactly zero. Other nodes do not affect this
 * element's base configuration and are deliberately excluded from the check.
 *
 * Component-wise equality avoids tolerance-based acceptance and underflow in
 * squared norms. NaN and infinite components also fail the zero comparison.
 * The field is read without modifying model geometry or material state.
 *
 * @param base_displacement Optional global nodal base translations and rotations.
 */
template<Index N>
void BeamElement<N>::validate_reference_base_displacement(const Field* base_displacement) const {
    // Check the complete local base before evaluating any linear operator.
    // An omitted field already represents zero displacement at every node.
    bool reference_base = true;
    if (base_displacement) {
        for (Index node : node_ids) {
            const Vec6 displacement = base_displacement->row_vec6(node);
            reference_base &= (displacement.array() == Precision(0)).all();
        }
    }
    logging::error(reference_base,
        "BeamElement: finite-rotation base displacement is not implemented yet for element ",
        this->label());
}

/**
 * Evaluates the currently supported beam operators through one mechanical entry.
 *
 * The existing beam formulations are linear and support only u0 = 0. The
 * ordinary elastic tangent is therefore evaluated at the reference state, while
 * the separate geometric stiffness is generated by the stress increment caused
 * by the supplied displacement u - u0 = u. Both an omitted base and an explicit
 * zero field represent u0 = 0. A nonzero base would require the missing
 * finite-rotation residual/tangent formulation and remains unsupported. Internal
 * forces use the same reference configuration through the linear law K u.
 */
template<Index N>
MapMatrix BeamElement<N>::evaluate(
    Precision*   tangent,
    Precision*   geometric_tangent,
    NodeData*    internal_force,
    const Field* target_displacement,
    const Field* target_temperature,
    const Field* base_displacement,
    const Field* base_temperature,
    bool         update_state
) {
    (void) update_state;

    validate_reference_base_displacement(base_displacement);
    logging::error(internal_force == nullptr || target_displacement != nullptr,
        "BeamElement: internal-force evaluation requires displacement for element ",
        this->label());

    if (tangent != nullptr) {
        MapMatrix mapped(tangent, N * 6, N * 6);
        mapped = stiffness_impl();
    }

    if (geometric_tangent != nullptr) {
        logging::error(target_displacement != nullptr,
            "BeamElement: geometric stiffness requires a target displacement from u0 = 0");
        MapMatrix mapped(geometric_tangent, N * 6, N * 6);
        mapped = stiffness_geom_impl(
            *target_displacement,
            target_temperature,
            base_temperature);
    }

    if (internal_force != nullptr) {
        // The linear beam uses
        //
        //     epsilon_mech = epsilon - epsilon_th
        //
        // for its axial constitutive response. In the linear element operator
        // this appears algebraically as the corresponding axial offset in K u.
        Eigen::Matrix<Precision, N * 6, 1> u_global;
        for (Index node = 0; node < N; ++node) {
            const Vec6 row =
                target_displacement->row_vec6(node_ids[node]);
            for (Index dof = 0; dof < 6; ++dof) {
                u_global(node * 6 + dof) = row(dof);
            }
        }

        StaticVector<N * 6> force = stiffness_impl() * u_global;

        auto material = get_material();
        if (target_temperature && material->has_thermal_expansion()) {
            const Precision T0    = material->get_thermal_zero_temperature();
            const Precision alpha = material->get_thermal_expansion();

            Precision temperature = Precision(0);
            for (Index node = 0; node < N; ++node) {
                temperature += (*target_temperature)(
                    node_ids[node], 0);
            }
            temperature /= static_cast<Precision>(N);

            const Precision axial_force =
                get_elasticity()->youngs * get_profile()->area_
                * alpha * (temperature - T0);

            StaticVector<N * 6> local_thermal = StaticVector<N * 6>::Zero();
            local_thermal(0)           = -axial_force;
            local_thermal((N - 1) * 6) =  axial_force;
            force -= transformation().transpose() * local_thermal;
        }

        for (Index node = 0; node < N; ++node) {
            const Index row = elem_nodal_offset + node;
            for (Index dof = 0; dof < 6; ++dof) {
                (*internal_force)(row, dof) += force(node * 6 + dof);
            }
        }
    }

    if (tangent != nullptr) {
        return MapMatrix(tangent, N * 6, N * 6);
    }
    if (geometric_tangent != nullptr) {
        return MapMatrix(geometric_tangent, N * 6, N * 6);
    }
    return MapMatrix(nullptr, 0, 0);
}

/**
 * Maps the formulation-specific beam mass matrix into caller-owned storage.
 *
 * @param buffer Caller-provided dense element-matrix storage.
 * @return Mapped beam mass matrix.
 */
template<Index N>
MapMatrix BeamElement<N>::mass(Precision* buffer) {
    MapMatrix result(buffer, N * 6, N * 6);
    result = mass_impl();
    return result;
}

template<Index N>
ElDofs BeamElement<N>::dofs() const {
    return ElDofs{true, true, true, true, true, true};
}

template<Index N>
Index BeamElement<N>::dimensions() const {
    return 3;
}

template<Index N>
Index BeamElement<N>::n_nodes() const {
    return N;
}

template<Index N>
Index BeamElement<N>::num_ip() const {
    return 1;
}

template<Index N>
const Index* BeamElement<N>::nodes() const {
    return node_ids.data();
}

template<Index N>
SurfacePtr BeamElement<N>::surface(Index surface_id) {
    (void) surface_id;
    return nullptr;
}

/**
 * Recovers linear beam section resultants in the original section frame.
 *
 * The base must be the implicit reference state or an explicit zero field at
 * the connected nodes. The target displacement supplies the total global beam
 * DOFs. The elastic force K u is corrected for axial thermal expansion, rotated
 * into the base section frame and written to element-nodal output rows using
 * the existing end-force sign convention. No base or material state is changed.
 *
 * @param section_forces Caller-owned element-nodal section force output.
 * @param target_displacement Total global nodal translations and rotations.
 * @param target_temperature Optional nodal temperatures for thermal expansion.
 * @param base_displacement Optional zero-valued reference displacement field.
 * @param base_temperature Unused by the reference-state linear beam law.
 * @return True after all element-nodal resultants have been written.
 */
template<Index N>
bool BeamElement<N>::compute_beam_section_forces(
    Field&       section_forces,
    const Field& target_displacement,
    const Field* target_temperature,
    const Field* base_displacement,
    const Field* base_temperature
) {
    // First element-nodal result row belonging to this element.
    Index offset = this->elem_nodal_offset;

    (void) base_temperature;
    validate_reference_base_displacement(base_displacement);
    Eigen::Matrix<Precision, N * 6, 1> u_global;
    for (Index i = 0; i < N; ++i) {
        const Index nid = node_ids[i];
        const Vec6 row = target_displacement.row_vec6(nid);
        for (Index d = 0; d < 6; ++d) {
            u_global(i * 6 + d) = row(d);
        }
    }

    const auto K_global          = stiffness_impl();
    const auto T_out             = transformation_base();
    StaticVector<N * 6> f_global = K_global * u_global;

    auto material = get_material();
    if (target_temperature && material->has_thermal_expansion()) {
        const Precision T0    = material->get_thermal_zero_temperature();
        const Precision alpha = material->get_thermal_expansion();
        Precision temperature = Precision(0);
        for (Index node = 0; node < N; ++node) {
            temperature += (*target_temperature)(
                node_ids[node], 0);
        }
        temperature /= static_cast<Precision>(N);

        const Precision axial_force =
            get_elasticity()->youngs * get_profile()->area_
            * alpha * (temperature - T0);

        StaticVector<N * 6> local_thermal = StaticVector<N * 6>::Zero();
        local_thermal(0)           = -axial_force;
        local_thermal((N - 1) * 6) =  axial_force;
        f_global -= transformation().transpose() * local_thermal;
    }

    const auto q_local = T_out * f_global;

    for (Index i = 0; i < N; ++i) {
        for (Index d = 0; d < 6; ++d) {
            section_forces(offset + i, d)
                = q_local(i * 6 + d) * ((i == 1 && N == 2) ? 1 : -1);
        }
    }
    return true;
}

template<Index N>
Precision BeamElement<N>::integrate_scalar_field(bool scale_by_density, const ScalarField& field) {
    const Precision L = length_reference();
    const Precision A = get_profile()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return Precision(0);

    // Sample the field on the reference axis so the integral is a dead load
    const Vec3 x_mid = midpoint_reference();

    Precision rho = Precision(1);
    if (scale_by_density) {
        auto mat = get_material();
        logging::error(mat && mat->has_density(),
            "BeamElement: material density is required when scale_by_density=true for element ", this->label());
        rho = mat->get_density();
    }
    return field(x_mid) * (rho * A * L);
}

template<Index N>
Vec3 BeamElement<N>::integrate_vector_field(bool scale_by_density, const VecField& field) {
    const Precision L = length_reference();
    const Precision A = get_profile()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return Vec3::Zero();

    // Sample the field on the reference axis so the integral is a dead load
    const Vec3 x_mid = midpoint_reference();

    Precision rho = Precision(1);
    if (scale_by_density) {
        auto mat = get_material();
        logging::error(mat && mat->has_density(),
            "BeamElement: material density is required when scale_by_density=true for element ", this->label());
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
    const Precision L = length_reference();
    const Precision A = get_profile()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return;

    // Sample the field on the reference axis so the integral is a dead load
    const Vec3 x_mid = midpoint_reference();

    Precision rho = Precision(1);
    if (scale_by_density) {
        auto mat = get_material();
        logging::error(mat && mat->has_density(),
            "BeamElement: material density is required when scale_by_density=true for element ", this->label());
        rho = mat->get_density();
    }

    const Vec3 F = field(x_mid) * (rho * A * L);
    const Precision share = Precision(1) / static_cast<Precision>(N);
    for (Index i = 0; i < N; ++i) {
        const Index n_id = node_ids[i];
        node_loads(n_id, 0) += share * F(0);
        node_loads(n_id, 1) += share * F(1);
        node_loads(n_id, 2) += share * F(2);
    }
}

/**
 * Integrates the scalar shape-function product operator of the field integration.
 *
 * The beam field integration samples the field at the reference nodal midpoint,
 * whose interpolation weights are 1/N, and distributes the resultant equally
 * to the N nodes. The consistent operator is therefore constant,
 *
 *     M(i,j) = rho A L0 / N^2.
 *
 * Degenerate beams without length or area contribute a zero operator, matching
 * the vanishing field integral.
 *
 * @param scale_by_density Multiply the measure by the material density.
 * @return Constant N x N operator in element node ordering.
 */
template<Index N>
DynamicMatrix BeamElement<N>::integrate_shape_products(bool scale_by_density) {
    const Precision L = length_reference();
    const Precision A = get_profile()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return DynamicMatrix::Zero(N, N);

    Precision rho = Precision(1);
    if (scale_by_density) {
        auto mat = get_material();
        logging::error(mat && mat->has_density(),
            "BeamElement: material density is required when scale_by_density=true for element ", this->label());
        rho = mat->get_density();
    }

    // Equal interpolation and distribution weights of the one-point rule
    const Precision share = Precision(1) / static_cast<Precision>(N);
    return DynamicMatrix::Constant(N, N, rho * A * L * share * share);
}

template<Index N>
Mat3 BeamElement<N>::integrate_tensor_field(bool scale_by_density, const TenField& field) {
    const Precision L = length_reference();
    const Precision A = get_profile()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return Mat3::Zero();

    // Sample the field on the reference axis so the integral is a dead load
    const Vec3 x_mid = midpoint_reference();

    Precision rho = Precision(1);
    if (scale_by_density) {
        auto mat = get_material();
        logging::error(mat && mat->has_density(),
            "BeamElement: material density is required when scale_by_density=true for element ", this->label());
        rho = mat->get_density();
    }
    return field(x_mid) * (rho * A * L);
}

} // namespace model
} // namespace fem
