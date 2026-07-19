/**
 * @file beam.inl
 * @brief Implements the common templated base for structural beam elements.
 *
 * This file contains the definitions declared by `BeamElement<N>` in
 * `beam.h`. It implements validated section and material access, beam geometry,
 * local section frames, rigid-axis offsets, matrix-interface adapters, result
 * recovery and midpoint field integration without introducing formulation-
 * specific stiffness or mass assumptions.
 *
 * The file is included by `beam.h` and must not be compiled as an independent
 * translation unit.
 *
 * @see beam.h
 */

namespace fem {
namespace model {

/**
 * @brief Constructs a beam element from its identifier and ordered connectivity.
 *
 * The constructor forwards the element identifier to `StructuralElement` and
 * stores both the fixed-size node connectivity and the optional orientation
 * node without further interpretation.
 */
template<Index N>
BeamElement<N>::BeamElement(
    ID                elem_id,
    std::array<ID, N> node_ids_in,
    ID                orientation_node_id
)
    : StructuralElement(elem_id)
    , node_ids(node_ids_in)
    , orientation_node_id_(orientation_node_id) {}

/**
 * @brief Destroys the beam element through the polymorphic base interface.
 */
template<Index N>
BeamElement<N>::~BeamElement() = default;

/**
 * @brief Returns the assigned section as a validated beam section.
 *
 * The function reports missing sections and incompatible section types through
 * the project logging facility before returning the typed section pointer.
 */
template<Index N>
BeamSection* BeamElement<N>::get_section() {
    logging::error(this->_section != nullptr,
        "Section not set for element ", this->elem_id);
    logging::error(this->_section->template as<BeamSection>() != nullptr,
        "Section is not a beam section for element ", this->elem_id);
    return this->_section->template as<BeamSection>();
}

/**
 * @brief Returns the geometric profile owned by the assigned beam section.
 */
template<Index N>
Profile* BeamElement<N>::get_profile() {
    return get_section()->profile_.get();
}

/**
 * @brief Returns the material assigned to the beam section.
 *
 * A missing material is reported before the section's shared material pointer
 * is returned.
 */
template<Index N>
material::MaterialPtr BeamElement<N>::get_material() {
    BeamSection* section = get_section();
    logging::error(section->material_ != nullptr,
        "Material not set for element ", this->elem_id);
    return section->material_;
}

/**
 * @brief Returns the isotropic elastic law assigned to the beam material.
 *
 * The complete section-material-elasticity chain is validated, including the
 * requirement that the active elasticity model is isotropic.
 */
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

/**
 * @brief Returns the optional global orientation-node identifier.
 */
template<Index N>
ID BeamElement<N>::orientation_node() const {
    return orientation_node_id_;
}

/**
 * @brief Reports whether an explicit orientation node is assigned.
 */
template<Index N>
bool BeamElement<N>::has_orientation_node() const {
    return orientation_node_id_ >= 0;
}

/**
 * @brief Determines and normalizes the section n1 orientation direction.
 *
 * An element-level orientation node takes precedence over a section-defined
 * direction. The function validates model-data access and rejects missing or
 * degenerate orientation definitions.
 */
template<Index N>
Vec3 BeamElement<N>::orientation_direction() {
    constexpr Precision kOrientationEps = static_cast<Precision>(1e-12);
    BeamSection* section = get_section();
    const bool section_has_direction = section && section->direction_.norm() > kOrientationEps;
    const bool element_has_node = has_orientation_node();

    logging::error(element_has_node || section_has_direction,
        "Beam element ", this->elem_id,
        " requires either an orientation node or a section-defined n1 vector");

    Vec3 n1_vec = Vec3::Zero();
    if (element_has_node) {
        logging::error(this->_model_data != nullptr,
            "no model data assigned to element ", this->elem_id);
        logging::error(this->_model_data->positions != nullptr,
            "positions field not set in model data");
        n1_vec = this->_model_data->positions->row_vec3(static_cast<Index>(orientation_node_id_));
    } else {
        n1_vec = section->direction_;
    }

    logging::error(n1_vec.norm() > kOrientationEps,
        "Orientation vector for element ", this->elem_id, " must be non-zero");
    return n1_vec.normalized();
}

/**
 * @brief Computes the polyline length through all consecutive beam nodes.
 */
template<Index N>
Precision BeamElement<N>::length() {
    Precision l = 0;
    for (Index i = 0; i < N - 1; i++) {
        const Vec3 x0 = this->node_position(static_cast<ID>(i));
        const Vec3 x1 = this->node_position(static_cast<ID>(i + 1));
        l += (x0 - x1).norm();
    }
    return l;
}

/**
 * @brief Computes the beam volume from section area and centerline length.
 */
template<Index N>
Precision BeamElement<N>::volume() {
    return get_profile()->area_ * length();
}

/**
 * @brief Builds the original local section frame from the beam axis and n1 direction.
 *
 * The local x-axis follows the first beam segment. The supplied orientation
 * direction initializes the local y-axis, while cross products orthogonalize
 * the final right-handed y-z basis.
 */
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

/**
 * @brief Computes the rotation from the base y-z frame to principal bending axes.
 *
 * The angle diagonalizes the section's two-dimensional inertia tensor. A
 * scale-aware tolerance suppresses numerically irrelevant product inertia.
 */
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

/**
 * @brief Builds the beam frame aligned with the principal bending axes.
 *
 * The base x-axis remains unchanged while the local y-z rows are rotated by
 * the principal-axis angle derived from the section inertia tensor.
 */
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

/**
 * @brief Returns the principal section frame used by internal beam matrices.
 */
template<Index N>
Mat3 BeamElement<N>::rotation_matrix() {
    return principal_rotation_matrix();
}

/**
 * @brief Builds the six-degree-of-freedom rigid-offset transformation.
 *
 * For a local point shift `d = (0, dy, dz)`, translational motion is augmented
 * by `theta x d` while the rotational components remain unchanged.
 */
template<Index N>
StaticMatrix<6, 6> BeamElement<N>::rigid_offset_6(Precision dy, Precision dz) {
    StaticMatrix<6, 6> B = StaticMatrix<6, 6>::Identity();

    // u += theta x d, with d=(0,dy,dz)
    // [ux] +=  ry*dz - rz*dy
    // [uy] += -rx*dz
    // [uz] +=  rx*dy
    B(0, 4) =  dz;
    B(0, 5) = -dy;

    B(1, 3) = -dz;
    B(2, 3) =  dy;

    return B;
}

/**
 * @brief Builds the block-diagonal rigid-offset transformation for all nodes.
 */
template<Index N>
StaticMatrix<N * 6, N * 6> BeamElement<N>::rigid_offset_N(Precision dy, Precision dz) {
    StaticMatrix<N * 6, N * 6> B = StaticMatrix<N * 6, N * 6>::Identity();
    const StaticMatrix<6, 6> Bi = rigid_offset_6(dy, dz);
    for (Index a = 0; a < N; ++a) {
        B.template block<6, 6>(a * 6, a * 6) = Bi;
    }
    return B;
}

/**
 * @brief Rotates a local y-z offset into the principal section frame.
 *
 * The operation applies the same in-plane rotation used to transform the base
 * section axes and writes the resulting components back through the references.
 */
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

/**
 * @brief Builds the block transformation from global to principal local DOFs.
 *
 * The same section rotation is inserted for the translational and rotational
 * triplets of every beam node.
 */
template<Index N>
StaticMatrix<N * 6, N * 6> BeamElement<N>::transformation() {
    StaticMatrix<N * 6, N * 6> T;
    T.setZero();
    Mat3 R = rotation_matrix();

    for (Index i = 0; i < N; i++) {
        for (Dim j = 0; j < 3; j++) {
            for (Dim k = 0; k < 3; k++) {
                T(i * 6 + j,     i * 6 + k)     = R(j, k);
                T(i * 6 + j + 3, i * 6 + k + 3) = R(j, k);
            }
        }
    }
    return T;
}

/**
 * @brief Builds the block transformation from global DOFs to the base section frame.
 *
 * Unlike `transformation()`, this variant does not apply the additional
 * principal-axis rotation and is therefore used for section-frame output.
 */
template<Index N>
StaticMatrix<N * 6, N * 6> BeamElement<N>::transformation_base() {
    StaticMatrix<N * 6, N * 6> T;
    T.setZero();
    Mat3 R = base_rotation_matrix();

    for (Index i = 0; i < N; i++) {
        for (Dim j = 0; j < 3; j++) {
            for (Dim k = 0; k < 3; k++) {
                T(i * 6 + j,     i * 6 + k)     = R(j, k);
                T(i * 6 + j + 3, i * 6 + k + 3) = R(j, k);
            }
        }
    }
    return T;
}

/**
 * @brief Writes the formulation-specific elastic stiffness into caller storage.
 */
template<Index N>
MapMatrix BeamElement<N>::stiffness(Precision* buffer) {
    MapMatrix result(buffer, N * 6, N * 6);
    result = stiffness_impl();
    return result;
}

/**
 * @brief Writes the formulation-specific geometric stiffness into caller storage.
 */
template<Index N>
MapMatrix BeamElement<N>::stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) {
    MapMatrix result(buffer, N * 6, N * 6);
    result = stiffness_geom_impl(ip_stress, ip_start_idx);
    return result;
}

/**
 * @brief Writes the formulation-specific mass matrix into caller storage.
 */
template<Index N>
MapMatrix BeamElement<N>::mass(Precision* buffer) {
    MapMatrix result(buffer, N * 6, N * 6);
    result = mass_impl();
    return result;
}

/**
 * @brief Reports that no generic nonlinear internal-force formulation exists.
 */
template<Index N>
void BeamElement<N>::compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) {
    (void)node_forces;
    (void)ip_stress;
    logging::error(false,
        "BeamElement: compute_internal_force_nonlinear is not implemented yet for element ", this->elem_id);
}

/**
 * @brief Recovers the concrete beam's stress state at its integration points.
 *
 * Empty integration-point sets are ignored. Otherwise the function delegates to
 * the formulation-specific stress/strain recovery with only the stress output.
 */
template<Index N>
void BeamElement<N>::compute_stress_state(
    Field&          stress_state,
    const Field&    displacement,
    int             offset,
    bool            use_green_lagrange_nl
) {
    RowMatrix rst = stress_strain_ip_rst();
    if (rst.rows() == 0) {
        return;
    }
    compute_stress_strain(nullptr, &stress_state, displacement, rst, offset, use_green_lagrange_nl);
}

/**
 * @brief Preserves the no-op base implementation for thermal loading.
 */
template<Index N>
void BeamElement<N>::apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) {
    (void)node_loads;
    (void)node_temp;
    (void)ref_temp;
}

/**
 * @brief Preserves the no-op base implementation for compliance evaluation.
 */
template<Index N>
void BeamElement<N>::compute_compliance(Field& displacement, Field& result) {
    (void)displacement;
    (void)result;
}

/**
 * @brief Preserves the no-op base implementation for compliance-angle derivatives.
 */
template<Index N>
void BeamElement<N>::compute_compliance_angle_derivative(Field& displacement, Field& result) {
    (void)displacement;
    (void)result;
}

/**
 * @brief Reports that generic beam shear-flow recovery is unavailable.
 */
template<Index N>
bool BeamElement<N>::compute_shear_flow(Field& shear_flow, const Field& displacement, int offset) {
    (void)shear_flow;
    (void)displacement;
    (void)offset;
    return false;
}

/**
 * @brief Returns the six active beam degrees of freedom per node.
 */
template<Index N>
ElDofs BeamElement<N>::dofs() const {
    return ElDofs{true, true, true, true, true, true};
}

/**
 * @brief Returns the three-dimensional embedding dimension.
 */
template<Index N>
Dim BeamElement<N>::dimensions() const {
    return 3;
}

/**
 * @brief Returns the compile-time number of beam nodes.
 */
template<Index N>
Dim BeamElement<N>::n_nodes() const {
    return N;
}

/**
 * @brief Returns the single integration point exposed by the base interface.
 */
template<Index N>
Dim BeamElement<N>::num_ip() const {
    return 1;
}

/**
 * @brief Returns contiguous access to the ordered beam connectivity.
 */
template<Index N>
const ID* BeamElement<N>::nodes() const {
    return node_ids.data();
}

/**
 * @brief Returns no surface representation for a generic beam element.
 */
template<Index N>
SurfacePtr BeamElement<N>::surface(ID surface_id) {
    (void)surface_id;
    return nullptr;
}

/**
 * @brief Recovers nodal beam resultants in the original section frame.
 *
 * Global nodal displacements are gathered in element order, multiplied by the
 * global element stiffness and transformed to the base section frame. The
 * established node-dependent sign convention is then applied to the result rows.
 */
template<Index N>
bool BeamElement<N>::compute_beam_section_forces(
    Field&       section_forces,
    const Field& displacement,
    int          offset
) {
    Eigen::Matrix<Precision, N * 6, 1> u_global;
    for (Index i = 0; i < N; ++i) {
        const ID nid = node_ids[i];
        // The displacement row follows the six-component beam ordering
        // [ux, uy, uz, rx, ry, rz].
        const Vec6 row = displacement.row_vec6(static_cast<Index>(nid));
        for (Index d = 0; d < 6; ++d) {
            u_global(i * 6 + d) = row(d);
        }
    }

    // global stiffness and local frame for output
    const auto K_global = stiffness_impl();
    // Output resultants are expressed in the base section frame.
    const auto T_out    = transformation_base();

    const auto f_global = K_global * u_global;
    // Keep the resultants associated with the element DOF reference line; no
    // additional shift to a section centroid is applied at this stage.
    const auto q_local  = T_out * f_global;

    for (Index i = 0; i < N; ++i) {

        for (Index d = 0; d < 6; ++d) {
            section_forces(static_cast<Index>(offset) + i, d)
              = q_local(i * 6 + d) * ((i == 1 && N == 2) ? 1 : -1);
        }
    }

    return true;
}

/**
 * @brief Reports that shell-section force recovery is not applicable to beams.
 */
template<Index N>
bool BeamElement<N>::compute_shell_section_forces(
    Field&       section_forces,
    Field&       contribution_count,
    const Field& displacement
) {
    (void)section_forces;
    (void)contribution_count;
    (void)displacement;
    return false;
}

/**
 * @brief Integrates a scalar field by midpoint sampling over the beam volume.
 *
 * The beam is represented by its physical length and profile area. The scalar
 * field is evaluated at the arithmetic mean of the beam-node positions, which
 * serves as the one-point sampling position for the current beam formulation.
 * The sampled value is multiplied by `A L` and, when requested, by the material
 * density. Degenerate beams or profiles contribute zero.
 *
 * @param scale_by_density Whether the material density is included in the
 *                         integration measure.
 * @param field Scalar field evaluated at the beam midpoint.
 * @return Volume integral of the scalar field, optionally density-scaled.
 */
template<Index N>
Precision BeamElement<N>::integrate_scalar_field(bool scale_by_density, const ScalarField& field) {
    // Compute the geometric measure used by the one-point beam integration.
    const Precision L = length();
    const Precision A = get_profile()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return Precision(0);

    // Use the arithmetic mean of all beam-node positions as the physical
    // sampling point for the constant one-point beam representation.
    Vec3 x_mid = Vec3::Zero();
    for (Index i = 0; i < N; ++i) x_mid += this->node_position(static_cast<ID>(i));
    x_mid /= static_cast<Precision>(N);

    // Determine the density factor. The default factor represents a purely
    // geometric volume integral; density scaling converts it to a mass-weighted
    // integral and requires a material density.
    Precision rho = 1.0;
    if (scale_by_density) {
        auto mat = get_material();
        logging::error(mat && mat->has_density(),
            "BeamElement: material density is required when scale_by_density=true for element ", this->elem_id);
        rho = mat->get_density();
    }

    // Apply the beam volume measure to the midpoint field value.
    return field(x_mid) * (rho * A * L);
}

/**
 * @brief Integrates a vector field by midpoint sampling over the beam volume.
 *
 * The vector field is evaluated once at the arithmetic mean of the beam-node
 * positions. Its value is multiplied by the beam volume `A L` and optionally
 * by the material density. This overload returns the total physical vector
 * integral and does not assemble nodal contributions.
 *
 * @param scale_by_density Whether the material density is included in the
 *                         integration measure.
 * @param field Vector field evaluated at the beam midpoint.
 * @return Volume integral of the vector field.
 */
template<Index N>
Vec3 BeamElement<N>::integrate_vector_field(bool scale_by_density, const VecField& field) {
    // Reject a beam with no physical length or no section area before sampling
    // the field or accessing a density factor.
    const Precision L = length();
    const Precision A = get_profile()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return Vec3::Zero();

    // Construct the physical midpoint used by the one-point integration rule.
    Vec3 x_mid = Vec3::Zero();
    for (Index i = 0; i < N; ++i) x_mid += this->node_position(static_cast<ID>(i));
    x_mid /= static_cast<Precision>(N);

    // Resolve the optional density weighting independently from the field
    // evaluation so the geometric and material measures remain explicit.
    Precision rho = 1.0;
    if (scale_by_density) {
        auto mat = get_material();
        logging::error(mat && mat->has_density(),
            "BeamElement: material density is required when scale_by_density=true for element ", this->elem_id);
        rho = mat->get_density();
    }

    // Return the midpoint value multiplied by the complete beam volume measure.
    return field(x_mid) * (rho * A * L);
}

/**
 * @brief Integrates a vector field and distributes the resultant to beam nodes.
 *
 * The field is sampled once at the arithmetic midpoint. The volume-integrated
 * vector is divided equally among all nodes and added to their translational
 * load components. Rotational DOFs are not populated by this generic load
 * integration path.
 *
 * @param node_loads Global nodal load field receiving the assembled force.
 * @param scale_by_density Whether the material density is included in the
 *                         integrated measure.
 * @param field Vector field evaluated at the beam midpoint.
 */
template<Index N>
void BeamElement<N>::integrate_vector_field(
    Field&          node_loads,
    bool            scale_by_density,
    const VecField& field
) {
    // Compute the physical beam volume and reject degenerate geometry before
    // modifying the caller-owned nodal load field.
    const Precision L = length();
    const Precision A = get_profile()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return;

    // Determine the single physical sampling point used by this midpoint rule.
    Vec3 x_mid = Vec3::Zero();
    for (Index i = 0; i < N; ++i) x_mid += this->node_position(static_cast<ID>(i));
    x_mid /= static_cast<Precision>(N);

    // Resolve the optional density factor before computing the total physical
    // force that will be distributed to the nodes.
    Precision rho = 1.0;
    if (scale_by_density) {
        auto mat = get_material();
        logging::error(mat && mat->has_density(),
            "BeamElement: material density is required when scale_by_density=true for element ", this->elem_id);
        rho = mat->get_density();
    }

    // Evaluate the total integrated force and assign an equal share to each
    // node's translational components in the global nodal load field.
    const Vec3 F = field(x_mid) * (rho * A * L);
    const Precision share = Precision(1) / static_cast<Precision>(N);
    for (Index i = 0; i < N; ++i) {
        const ID n_id = node_ids[i];

        // Assemble only translational components; no generic moment resultant
        // is inferred from the distributed vector field.
        node_loads(n_id, 0) += share * F(0);
        node_loads(n_id, 1) += share * F(1);
        node_loads(n_id, 2) += share * F(2);
    }
}

/**
 * @brief Integrates a second-order tensor field by midpoint sampling.
 *
 * The tensor field is evaluated once at the arithmetic beam midpoint and
 * multiplied by the physical beam volume `A L`. Optional density scaling uses
 * the material density after validating that it is available. This overload
 * returns the total tensor integral and does not assemble nodal values.
 *
 * @param scale_by_density Whether the material density is included in the
 *                         integration measure.
 * @param field Tensor field evaluated at the beam midpoint.
 * @return Volume integral of the tensor field.
 */
template<Index N>
Mat3 BeamElement<N>::integrate_tensor_field(bool scale_by_density, const TenField& field) {
    // Reject a degenerate beam before forming its physical integration measure.
    const Precision L = length();
    const Precision A = get_profile()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return Mat3::Zero();

    // Gather the nodal positions once and construct the midpoint used by the
    // one-point tensor integration rule.
    Vec3 x_mid = Vec3::Zero();
    for (Index i = 0; i < N; ++i) x_mid += this->node_position(static_cast<ID>(i));
    x_mid /= static_cast<Precision>(N);

    // Resolve the optional density weighting while keeping the geometric
    // measure explicit in the final tensor integral.
    Precision rho = 1.0;
    if (scale_by_density) {
        auto mat = get_material();
        logging::error(mat && mat->has_density(),
            "BeamElement: material density is required when scale_by_density=true for element ", this->elem_id);
        rho = mat->get_density();
    }

    // Apply the complete density-scaled beam volume to the sampled tensor.
    return field(x_mid) * (rho * A * L);
}

} // namespace model
} // namespace fem
