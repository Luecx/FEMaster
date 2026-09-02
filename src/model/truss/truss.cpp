/**
 * @file truss.cpp
 * @brief Implements the three-dimensional two-node truss element.
 *
 * The truss evaluates axial material response at one globally enumerated
 * material point. Linearized output uses Cauchy stress; the nonlinear element
 * uses Green-Lagrange strain and PK2 stress in a Total-Lagrangian formulation.
 * Nonlinear internal force and tangent are assembled from one constitutive trial
 * evaluation. State-neutral operators read committed history without storing a
 * constitutive output state.
 *
 * Constitutive tangents follow the nullable-pointer `Elasticity::evaluate`
 * contract. Stiffness operators request a tangent explicitly, while prestress
 * recovery and result output request stress only. Residual-only nonlinear
 * assembly therefore propagates a null tangent pointer into the material model.
 *
 * @see T3
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#include "truss.h"

#include <cmath>

namespace fem {
namespace model {
namespace {

/**
 * Returns the geometric midpoint of the truss in the current configuration.
 *
 * The midpoint is used as the single sampling location for distributed scalar,
 * vector and tensor fields integrated over the current truss volume.
 *
 * @param elem Truss whose current midpoint is requested.
 * @return Current geometric midpoint.
 */
Vec3 midpoint(T3& elem) {
    return (elem.node_position_current(0) + elem.node_position_current(1)) * Precision(0.5);
}

/**
 * Returns the multiplicative density factor for field integration.
 *
 * Unscaled integration uses unity. Density-scaled integration requires a
 * material density and returns that value so the calling integration routine can
 * convert volume-based quantities to mass-based quantities.
 *
 * @param elem Truss providing the material definition.
 * @param scale_by_density Select density-scaled or purely geometric integration.
 * @return Unity or the assigned material density.
 */
Precision density_scale(T3& elem, bool scale_by_density) {
    if (!scale_by_density) {
        return Precision(1);
    }

    auto mat = elem.get_material();
    logging::error(mat && mat->has_density(),
        "T3: material density is required when scale_by_density=true for element ", elem.elem_id);
    return mat->get_density();
}

} // namespace

/**
 * Constructs a two-node truss from its element id and global node ids.
 *
 * Section, material, model data and material-point storage are bound later by
 * the compiled model infrastructure.
 *
 * @param elem_id Global element identifier.
 * @param node_ids_in Global identifiers of the two truss nodes.
 */
T3::T3(ID elem_id, std::array<ID, N> node_ids_in)
    : StructuralElement(elem_id),
      node_ids(node_ids_in) {}

ElDofs T3::dofs() const {
    return ElDofs{true, true, true, false, false, false};
}

Dim T3::dimensions() const {
    return 3;
}

Dim T3::n_nodes() const {
    return N;
}

Dim T3::num_ip() const {
    return 1;
}

const ID* T3::nodes() const {
    return node_ids.data();
}

SurfacePtr T3::surface(ID surface_id) {
    (void) surface_id;
    return nullptr;
}

std::string T3::type_name() const {
    return "T3";
}

/**
 * Resolves and validates the truss section assigned to the element.
 *
 * @return Assigned `TrussSection`.
 */
TrussSection* T3::get_section() {
    logging::error(this->_section != nullptr,
        "T3: missing section for element ", this->elem_id);

    auto* section = this->_section->template as<TrussSection>();
    logging::error(section != nullptr,
        "T3: section is not a truss section for element ", this->elem_id);
    return section;
}

/**
 * Resolves the material referenced by the assigned truss section.
 *
 * @return Material assigned through the truss section.
 */
material::MaterialPtr T3::get_material() {
    TrussSection* section = get_section();
    logging::error(section->material_ != nullptr,
        "T3: no material set for element ", this->elem_id);
    return section->material_;
}

/**
 * Resolves the elastic constitutive law required by the truss formulation.
 *
 * @return Non-owning pointer to the assigned elasticity model.
 */
material::Elasticity* T3::get_elasticity() {
    auto mat = get_material();
    logging::error(mat->has_elasticity(),
        "T3: material has no elasticity for element ", this->elem_id);
    return mat->elasticity().get();
}

/**
 * Returns one nodal position from the undeformed reference configuration.
 *
 * @param local_node Local node index zero or one.
 * @return Reference position of the requested node.
 */
Vec3 T3::node_position_reference(Index local_node) const {
    logging::error(local_node < N,
        "T3: local node index out of range in element ", this->elem_id);
    logging::error(this->_model_data != nullptr,
        "T3: no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions_reference != nullptr,
        "T3: reference positions field not set in model data");

    return this->_model_data->positions_reference->row_vec3(static_cast<Index>(node_ids[local_node]));
}

/**
 * Returns one nodal position from the current configuration.
 *
 * @param local_node Local node index zero or one.
 * @return Current position of the requested node.
 */
Vec3 T3::node_position_current(Index local_node) const {
    logging::error(local_node < N,
        "T3: local node index out of range in element ", this->elem_id);
    logging::error(this->_model_data != nullptr,
        "T3: no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions != nullptr,
        "T3: current positions field not set in model data");

    return this->_model_data->positions->row_vec3(static_cast<Index>(node_ids[local_node]));
}

Precision T3::length_reference() const {
    return (node_position_reference(1) - node_position_reference(0)).norm();
}

Precision T3::length_current() const {
    return (node_position_current(1) - node_position_current(0)).norm();
}

/**
 * Returns the axial stretch of the current configuration.
 *
 *     lambda = l / L0.
 *
 * @return Positive axial stretch.
 */
Precision T3::stretch() const {
    const Precision L0 = length_reference();
    const Precision l  = length_current();

    logging::error(L0 > Precision(0),
        "T3: zero reference length for element ", this->elem_id);
    logging::error(l > Precision(0),
        "T3: zero current length for element ", this->elem_id);

    return l / L0;
}

/**
 * Returns the unit vector along the undeformed truss axis.
 *
 * @return Reference axial direction.
 */
Vec3 T3::direction_reference() const {
    const Precision L0 = length_reference();
    logging::error(L0 > Precision(0),
        "T3: zero reference length in element ", this->elem_id);

    return (node_position_reference(1) - node_position_reference(0)) / L0;
}

/**
 * Returns the unit vector along the current truss axis.
 *
 * @return Current axial direction.
 */
Vec3 T3::direction_current() const {
    const Precision l = length_current();
    logging::error(l > Precision(0),
        "T3: zero current length in element ", this->elem_id);

    return (node_position_current(1) - node_position_current(0)) / l;
}

Precision T3::length() {
    return length_current();
}

Vec3 T3::direction() {
    return direction_current();
}

/**
 * Returns the current geometric volume `A0 l` of the truss.
 *
 * @return Current line volume based on reference area.
 */
Precision T3::volume() {
    return get_section()->area_ * length_current();
}

/**
 * Assembles the linear axial stiffness in the reference configuration.
 *
 * For the linearized axial material tangent `C`, reference area `A0`, reference
 * length `L0` and reference direction `n0`, the translational block is
 *
 *     k = A0 C / L0 (n0 n0^T).
 *
 * The material tangent is required by this operator and is therefore requested
 * explicitly through a non-null tangent pointer. Constitutive state remains
 * state-neutral because no output state is supplied.
 *
 * @param buffer Caller-provided six-by-six element-matrix storage.
 * @return Mapped linear stiffness matrix.
 */
MapMatrix T3::stiffness(Precision* buffer) {
    // Collect reference geometry used by the complete linear operator.
    const Precision A0 = get_section()->area_;
    const Precision L0 = length_reference();
    const Vec3      n0 = direction_reference();

    logging::error(L0 > Precision(0),
        "T3: zero reference length in stiffness for element ", this->elem_id);

    // Evaluate the material tangent at zero infinitesimal strain. Stress is an
    // unavoidable output of the common constitutive interface but is not used by
    // this purely material stiffness operator.
    const AxialStrainLinearized axial_strain(Precision(0));
    AxialStressCauchy           axial_stress;
    Precision                   material_tangent = Precision(0);

    auto elasticity = get_elasticity();
    logging::error(elasticity->supports_axial_linearized(),
        "T3: material does not support linearized axial evaluation for element ", this->elem_id);

    const Index      state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);

    elasticity->evaluate(
        axial_strain,
        old_state,
        nullptr,
        axial_stress,
        &material_tangent
    );

    // Embed the scalar axial stiffness into the three translational directions.
    const Mat3 k = (A0 * material_tangent / L0) * (n0 * n0.transpose());

    StaticMatrix<N * 3, N * 3> K = StaticMatrix<N * 3, N * 3>::Zero();
    K.block(0, 0, 3, 3) =  k;
    K.block(0, 3, 3, 3) = -k;
    K.block(3, 0, 3, 3) = -k;
    K.block(3, 3, 3, 3) =  k;

    MapMatrix result(buffer, N * 3, N * 3);
    result = K;
    return result;
}

/**
 * Assembles the geometric stiffness associated with a supplied displacement
 * state without modifying persistent material history.
 *
 * The displacement difference is projected onto the reference truss axis to
 * obtain infinitesimal axial strain. Only the resulting Cauchy stress is needed:
 *
 *     K_geo = A0 sigma / L0 I.
 *
 * The material tangent is therefore deliberately omitted from this constitutive
 * call.
 *
 * @param buffer Caller-provided six-by-six element-matrix storage.
 * @param displacement Global nodal displacement field defining the prestress.
 * @return Mapped geometric stiffness matrix.
 */
MapMatrix T3::stiffness_geom(Precision* buffer, const Field& displacement) {
    const Precision L0 = length_reference();
    logging::error(L0 > Precision(0),
        "T3: zero reference length in stiffness_geom for element ", this->elem_id);

    // Reconstruct axial prestress from the supplied displacement state.
    const Vec3 u0 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
    const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[1]));

    const AxialStrainLinearized axial_strain((u1 - u0).dot(direction_reference()) / L0);
    AxialStressCauchy           axial_stress;

    auto elasticity = get_elasticity();
    logging::error(elasticity->supports_axial_linearized(),
        "T3: material does not support linearized axial evaluation for element ", this->elem_id);

    // Prestress assembly needs sigma but not d sigma/d epsilon.
    const Index      state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
    elasticity->evaluate(axial_strain, old_state, nullptr, axial_stress, nullptr);

    // The stress contribution multiplies the three-dimensional identity block.
    const Mat3 k = (get_section()->area_ * axial_stress.value() / L0) * Mat3::Identity();

    StaticMatrix<N * 3, N * 3> K = StaticMatrix<N * 3, N * 3>::Zero();
    K.block(0, 0, 3, 3) =  k;
    K.block(0, 3, 3, 3) = -k;
    K.block(3, 0, 3, 3) = -k;
    K.block(3, 3, 3, 3) =  k;

    MapMatrix result(buffer, N * 3, N * 3);
    result = K;
    return result;
}

/**
 * Evaluates nonlinear Total-Lagrangian equilibrium and optionally assembles the
 * consistent truss tangent from one axial constitutive update.
 *
 * Green-Lagrange strain follows from the current stretch. The material reads the
 * committed state and writes the physical trial state. The PK2 stress produces
 *
 *     f_int = A0 lambda S n.
 *
 * If matrix storage is supplied, the same constitutive candidate additionally
 * provides the consistent scalar derivative and the tangent blocks
 *
 *     K_mat = A0 C_alg lambda^2 / L0 (n n^T)
 *     K_geo = A0 S / L0 I.
 *
 * For `buffer == nullptr`, the material receives `tangent == nullptr`; stress and
 * trial state are still updated identically, but no constitutive tangent is
 * constructed.
 *
 * @param buffer Optional six-by-six tangent storage; null requests force only.
 * @param nodal_forces Global nodal internal-force field to increment.
 * @param displacement Trial displacement field; current positions already carry
 *                     the configuration used by the truss kinematics.
 * @return Mapped tangent matrix, or an empty map for residual-only evaluation.
 */
MapMatrix T3::stiffness_tangent(Precision*   buffer,
                                NodeData&    nodal_forces,
                                const Field& displacement) {
    (void) displacement;

    logging::error(nodal_forces.components >= 3,
        "T3: nonlinear internal force requires at least three nodal components");

    // Evaluate current kinematics in the Total-Lagrangian formulation.
    const Precision A0     = get_section()->area_;
    const Precision L0     = length_reference();
    const Precision lambda = stretch();

    logging::error(L0 > Precision(0),
        "T3: zero reference length in nonlinear tangent for element ", this->elem_id);

    const AxialStrainGreenLagrange strain =
        AxialStrainGreenLagrange::from_stretch(lambda);
    AxialStressPK2 stress;
    Precision      material_tangent = Precision(0);

    auto elasticity = get_elasticity();
    logging::error(elasticity->supports_axial_green_lagrange(),
        "T3: material does not support Green-Lagrange axial evaluation for element ", this->elem_id);

    // Perform one committed -> trial constitutive update. Tangent construction is
    // tied directly to whether the element matrix itself is requested.
    const Index      state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
    Precision*       new_state = &(*this->_model_data->material_state_new)(state_row, 0);

    elasticity->evaluate(
        strain,
        old_state,
        new_state,
        stress,
        buffer != nullptr ? &material_tangent : nullptr
    );

    // Scatter equal and opposite current-direction internal forces.
    const Vec3  n      = direction_current();
    const Vec3  force  = A0 * lambda * stress.value() * n;
    const Index node_0 = static_cast<Index>(node_ids[0]);
    const Index node_1 = static_cast<Index>(node_ids[1]);

    for (Dim d = 0; d < 3; ++d) {
        nodal_forces(node_0, d) -= force(d);
        nodal_forces(node_1, d) += force(d);
    }

    // Residual-only evaluation ends after the complete physical material update.
    if (buffer == nullptr) {
        return MapMatrix(nullptr, 0, 0);
    }

    // Assemble material and geometric contributions from the same stress/tangent
    // pair that produced the internal force.
    const Mat3 material_block =
        (A0 * material_tangent * lambda * lambda / L0) * (n * n.transpose());
    const Mat3 geometric_block =
        (A0 * stress.value() / L0) * Mat3::Identity();
    const Mat3 block = material_block + geometric_block;

    StaticMatrix<N * 3, N * 3> tangent = StaticMatrix<N * 3, N * 3>::Zero();
    tangent.block(0, 0, 3, 3) =  block;
    tangent.block(0, 3, 3, 3) = -block;
    tangent.block(3, 0, 3, 3) = -block;
    tangent.block(3, 3, 3, 3) =  block;

    MapMatrix mapped(buffer, N * 3, N * 3);
    mapped = tangent;
    return mapped;
}

/**
 * Assembles the lumped translational mass matrix of the truss.
 *
 * For material density `rho`, reference area `A0` and reference length `L0`, the
 * total mass is
 *
 *     m = rho A0 L0.
 *
 * Half is assigned to each node and repeated on its three translational DOFs.
 * Materials without density produce a zero matrix.
 *
 * @param buffer Caller-provided matrix storage.
 * @return Mapped six-by-six lumped mass matrix.
 */
MapMatrix T3::mass(Precision* buffer) {
    StaticMatrix<N * 3, N * 3> M = StaticMatrix<N * 3, N * 3>::Zero();

    auto mat = get_material();
    if (mat->has_density()) {
        const Precision rho = mat->get_density();
        const Precision A   = get_section()->area_;
        const Precision L0  = length_reference();
        const Precision m   = rho * A * L0;

        for (Index i = 0; i < N; ++i) {
            M.block(i * 3, i * 3, 3, 3) = Mat3::Identity() * (m / Precision(2));
        }
    }

    MapMatrix result(buffer, N * 3, N * 3);
    result = M;
    return result;
}

RowMatrix T3::stress_strain_nodal_rst() {
    RowMatrix rst(N, 3);
    rst.setZero();
    rst(0, 0) = Precision(-1);
    rst(1, 0) = Precision(1);
    return rst;
}

RowMatrix T3::stress_strain_ip_rst() {
    RowMatrix rst(1, 3);
    rst.setZero();
    return rst;
}

/**
 * Integrates a scalar field over the current truss volume using midpoint
 * evaluation.
 *
 * @param scale_by_density Multiply by material density when true.
 * @param field Spatial scalar field.
 * @return Integrated scalar quantity.
 */
Precision T3::integrate_scalar_field(bool scale_by_density, const ScalarField& field) {
    const Precision L = length_current();
    const Precision A = get_section()->area_;

    if (L <= Precision(0) || A <= Precision(0)) {
        return Precision(0);
    }

    return field(midpoint(*this)) * density_scale(*this, scale_by_density) * A * L;
}

/**
 * Integrates a vector field over the current truss volume using midpoint
 * evaluation.
 *
 * @param scale_by_density Multiply by material density when true.
 * @param field Spatial vector field.
 * @return Integrated vector quantity.
 */
Vec3 T3::integrate_vector_field(bool scale_by_density, const VecField& field) {
    const Precision L = length_current();
    const Precision A = get_section()->area_;

    if (L <= Precision(0) || A <= Precision(0)) {
        return Vec3::Zero();
    }

    return field(midpoint(*this)) * density_scale(*this, scale_by_density) * A * L;
}

/**
 * Integrates a distributed vector field and scatters the equivalent force equally
 * to both truss nodes.
 *
 * @param node_loads Global nodal load field to increment.
 * @param scale_by_density Multiply by material density when true.
 * @param field Spatial vector field.
 */
void T3::integrate_vector_field(Field& node_loads, bool scale_by_density, const VecField& field) {
    const Precision L = length_current();
    const Precision A = get_section()->area_;

    if (L <= Precision(0) || A <= Precision(0)) {
        return;
    }

    const Vec3 force = field(midpoint(*this)) * density_scale(*this, scale_by_density) * A * L;

    for (Index i = 0; i < N; ++i) {
        const Index node = static_cast<Index>(node_ids[i]);
        node_loads(node, 0) += force(0) * Precision(0.5);
        node_loads(node, 1) += force(1) * Precision(0.5);
        node_loads(node, 2) += force(2) * Precision(0.5);
    }
}

/**
 * Integrates a second-order tensor field over the current truss volume.
 *
 * @param scale_by_density Multiply by material density when true.
 * @param field Spatial tensor field.
 * @return Integrated tensor quantity.
 */
Mat3 T3::integrate_tensor_field(bool scale_by_density, const TenField& field) {
    const Precision L = length_current();
    const Precision A = get_section()->area_;

    if (L <= Precision(0) || A <= Precision(0)) {
        return Mat3::Zero();
    }

    return field(midpoint(*this)) * density_scale(*this, scale_by_density) * A * L;
}

/**
 * Applies equivalent nodal loading from a prescribed temperature field.
 *
 * Thermal expansion is currently not implemented for T3. The function therefore
 * remains a deliberate no-op while satisfying the structural-element interface.
 */
void T3::apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) {
    (void) node_loads;
    (void) node_temp;
    (void) ref_temp;
}

/**
 * Recovers axial strain and Cauchy stress at requested truss output locations.
 *
 * Linearized recovery projects relative displacement onto the reference axis.
 * Nonlinear recovery derives Green-Lagrange strain from stretch, evaluates PK2
 * stress and converts it to physical Cauchy stress through
 *
 *     sigma = lambda S.
 *
 * Result recovery never uses a constitutive tangent and therefore passes
 * `nullptr` explicitly for tangent output in both branches.
 *
 * @param strain Optional output field receiving axial strain in component zero.
 * @param stress Optional output field receiving axial Cauchy stress in component zero.
 * @param displacement Global nodal displacement field for linearized recovery.
 * @param rst Requested natural output coordinates.
 * @param offset First output row belonging to this element.
 * @param use_green_lagrange_nl Select finite-strain or linearized recovery.
 */
void T3::compute_stress_strain(Field*           strain,
                               Field*           stress,
                               const Field&     displacement,
                               const RowMatrix& rst,
                               int              offset,
                               bool             use_green_lagrange_nl) {
    logging::error(strain != nullptr || stress != nullptr,
        "T3: compute_stress_strain requires at least one output field");
    logging::error(rst.cols() >= 1,
        "T3: stress/strain coordinates require at least 1 column");

    Precision strain_value = Precision(0);
    Precision stress_value = Precision(0);

    const Index      state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);

    if (use_green_lagrange_nl) {
        // Recover the finite-strain work-conjugate material pair and convert PK2
        // stress to physical Cauchy stress for output.
        const Precision lambda = stretch();
        const AxialStrainGreenLagrange axial_strain =
            AxialStrainGreenLagrange::from_stretch(lambda);
        AxialStressPK2 axial_stress;

        auto elasticity = get_elasticity();
        logging::error(elasticity->supports_axial_green_lagrange(),
            "T3: material does not support Green-Lagrange axial evaluation for element ", this->elem_id);

        elasticity->evaluate(axial_strain, old_state, nullptr, axial_stress, nullptr);

        strain_value = axial_strain.value();
        stress_value = lambda * axial_stress.value();
    } else {
        // Recover infinitesimal axial strain from the reference-axis displacement
        // difference and evaluate Cauchy stress directly.
        const Precision L0 = length_reference();
        logging::error(L0 > Precision(0),
            "T3: zero reference length in compute_stress_strain for element ", this->elem_id);

        const Vec3 u0 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
        const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[1]));

        const AxialStrainLinearized axial_strain((u1 - u0).dot(direction_reference()) / L0);
        AxialStressCauchy           axial_stress;

        auto elasticity = get_elasticity();
        logging::error(elasticity->supports_axial_linearized(),
            "T3: material does not support linearized axial evaluation for element ", this->elem_id);

        elasticity->evaluate(axial_strain, old_state, nullptr, axial_stress, nullptr);

        strain_value = axial_strain.value();
        stress_value = axial_stress.value();
    }

    // Replicate the constant axial state to every requested output coordinate.
    for (Index i = 0; i < static_cast<Index>(rst.rows()); ++i) {
        const Index row = static_cast<Index>(offset) + i;

        if (strain) {
            for (Index j = 0; j < strain->components; ++j) {
                (*strain)(row, j) = Precision(0);
            }
            (*strain)(row, 0) = strain_value;
        }

        if (stress) {
            for (Index j = 0; j < stress->components; ++j) {
                (*stress)(row, j) = Precision(0);
            }
            (*stress)(row, 0) = stress_value;
        }
    }
}

/**
 * Computes the element compliance contribution `u^T K u`.
 *
 * @param displacement Global nodal displacement field.
 * @param result Element result field receiving the scalar compliance.
 */
void T3::compute_compliance(Field& displacement, Field& result) {
    Precision buffer[N * 3 * N * 3] {};
    MapMatrix K = stiffness(buffer);

    StaticVector<N * 3> u;
    for (Index i = 0; i < N; ++i) {
        const Vec3 ui = displacement.row_vec3(static_cast<Index>(node_ids[i]));
        for (Index d = 0; d < 3; ++d) {
            u(i * 3 + d) = ui(d);
        }
    }

    result(static_cast<Index>(this->elem_id), 0) = u.dot(K * u);
}

/**
 * Recovers the linearized axial section force for beam-style result output.
 *
 * Only Cauchy stress is required. The constitutive tangent is therefore omitted,
 * and axial force follows directly from
 *
 *     N = A0 sigma.
 *
 * @param section_forces Element-nodal section-force output field.
 * @param displacement Global nodal displacement field.
 * @param offset First element-nodal output row belonging to this truss.
 * @return Always true after both nodal result rows were written.
 */
bool T3::compute_beam_section_forces(Field&       section_forces,
                                     const Field& displacement,
                                     int          offset) {
    // Reconstruct infinitesimal axial strain from the reference geometry.
    const Vec3 u0 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
    const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[1]));

    const Precision L0 = length_reference();
    logging::error(L0 > Precision(0),
        "T3: zero reference length in compute_beam_section_forces for element ", this->elem_id);

    const AxialStrainLinearized axial_strain((u1 - u0).dot(direction_reference()) / L0);
    AxialStressCauchy           axial_stress;

    auto elasticity = get_elasticity();
    logging::error(elasticity->supports_axial_linearized(),
        "T3: material does not support linearized axial evaluation for element ", this->elem_id);

    // Section-force recovery requires stress only and remains state-neutral.
    const Index      state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
    elasticity->evaluate(axial_strain, old_state, nullptr, axial_stress, nullptr);

    const Precision axial_force = get_section()->area_ * axial_stress.value();

    // The T3 has a constant axial force; copy it to both element-nodal rows and
    // clear all unsupported section-force components.
    for (Index i = 0; i < N; ++i) {
        for (Index d = 0; d < section_forces.components; ++d) {
            section_forces(static_cast<Index>(offset) + i, d) = Precision(0);
        }
        section_forces(static_cast<Index>(offset) + i, 0) = axial_force;
    }

    return true;
}

} // namespace model
} // namespace fem