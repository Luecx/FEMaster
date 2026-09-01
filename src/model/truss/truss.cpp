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

Vec3 midpoint(T3& elem) {
    return (elem.node_position_current(0) + elem.node_position_current(1)) * Precision(0.5);
}

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

TrussSection* T3::get_section() {
    logging::error(this->_section != nullptr,
        "T3: missing section for element ", this->elem_id);

    auto* section = this->_section->template as<TrussSection>();
    logging::error(section != nullptr,
        "T3: section is not a truss section for element ", this->elem_id);
    return section;
}

material::MaterialPtr T3::get_material() {
    TrussSection* section = get_section();
    logging::error(section->material_ != nullptr,
        "T3: no material set for element ", this->elem_id);
    return section->material_;
}

material::Elasticity* T3::get_elasticity() {
    auto mat = get_material();
    logging::error(mat->has_elasticity(),
        "T3: material has no elasticity for element ", this->elem_id);
    return mat->elasticity().get();
}

Vec3 T3::node_position_reference(Index local_node) const {
    logging::error(local_node < N,
        "T3: local node index out of range in element ", this->elem_id);
    logging::error(this->_model_data != nullptr,
        "T3: no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions_reference != nullptr,
        "T3: reference positions field not set in model data");

    return this->_model_data->positions_reference->row_vec3(static_cast<Index>(node_ids[local_node]));
}

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

Precision T3::stretch() const {
    const Precision L0 = length_reference();
    const Precision l  = length_current();

    logging::error(L0 > Precision(0),
        "T3: zero reference length for element ", this->elem_id);
    logging::error(l > Precision(0),
        "T3: zero current length for element ", this->elem_id);

    return l / L0;
}

Vec3 T3::direction_reference() const {
    const Precision L0 = length_reference();
    logging::error(L0 > Precision(0),
        "T3: zero reference length in element ", this->elem_id);

    return (node_position_reference(1) - node_position_reference(0)) / L0;
}

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

Precision T3::volume() {
    return get_section()->area_ * length_current();
}

/**
 * Assembles the linear axial stiffness in the reference configuration.
 *
 * The linear truss operator is based on infinitesimal axial strain, Cauchy
 * stress and the undeformed element axis. For an axial material tangent `C`,
 * the three-dimensional nodal block is
 *
 *     k = A0 C/L0 (n0 tensor n0),
 *
 * where `A0`, `L0` and `n0` are the reference area, length and unit direction.
 * The usual two-node difference operator expands this block into the complete
 * six-by-six element matrix.
 *
 * This function deliberately contains no finite-strain kinematics. Current
 * stretch, current direction, Green-Lagrange strain and PK2 stress belong to
 * `stiffness_tangent()`, which is the nonlinear equilibrium operator.
 *
 * The constitutive call is state-neutral: committed history may be inspected,
 * but no target state is supplied and persistent trial history is not modified.
 *
 * @param buffer Caller-provided six-by-six element-matrix storage.
 * @return Mapped linear stiffness in global translational DOF ordering.
 */
MapMatrix T3::stiffness(Precision* buffer) {
    // Linear stiffness is defined entirely in the undeformed reference
    // configuration. The area, length and direction therefore all refer to the
    // reference truss geometry and do not depend on the current nodal positions.
    const Precision A0 = get_section()->area_;
    const Precision L0 = length_reference();
    const Vec3      n0 = direction_reference();

    logging::error(L0 > Precision(0),
        "T3: zero reference length in stiffness for element ", this->elem_id);

    // Query the linearized axial constitutive tangent. Zero strain is sufficient
    // because stiffness() represents the linear material operator itself rather
    // than the tangent of a finite-deformation trial configuration.
    const AxialStrainLinearized axial_strain(Precision(0));
    AxialStressCauchy axial_stress;
    Precision         material_tangent = Precision(0);

    auto elasticity = get_elasticity();
    logging::error(elasticity->supports_axial_linearized(),
        "T3: material does not support linearized axial evaluation for element ", this->elem_id);

    // The linear operator may read committed material history, but it must not
    // create or advance a trial state. A null target state makes this read-only
    // constitutive intent explicit at the call site.
    const Index      state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
    elasticity->evaluate(axial_strain, old_state, nullptr, axial_stress, material_tangent);

    // Embed the scalar axial stiffness A0*C/L0 into the three translational DOFs
    // using the dyadic product of the reference unit direction with itself.
    const Mat3 k = (A0 * material_tangent / L0) * (n0 * n0.transpose());

    // Apply the two-node displacement difference operator. Equal diagonal blocks
    // and opposite off-diagonal blocks produce equal and opposite nodal forces
    // while retaining only deformation along the reference truss axis.
    StaticMatrix<N * 3, N * 3> K = StaticMatrix<N * 3, N * 3>::Zero();
    K.block(0, 0, 3, 3) =  k;
    K.block(0, 3, 3, 3) = -k;
    K.block(3, 0, 3, 3) = -k;
    K.block(3, 3, 3, 3) =  k;

    // Copy the fixed-size local matrix into caller-owned storage and expose it
    // through the common dynamically mapped element-matrix interface.
    MapMatrix result(buffer, N * 3, N * 3);
    result = K;
    return result;
}

/**
 * Assembles the geometric stiffness associated with a supplied displacement
 * state without modifying persistent material history.
 *
 * The displacement difference is projected onto the reference truss axis to
 * obtain the infinitesimal axial strain used for prestress evaluation. The
 * corresponding axial Cauchy stress supplies the stress-dependent block
 *
 *     K_geo = A0 sigma/L0 I.
 *
 * The constitutive evaluation reads committed history and receives no target
 * state. The displacement field is therefore the complete external input needed
 * for prestress stiffness assembly.
 *
 * @param buffer Caller-provided six-by-six element-matrix storage.
 * @param displacement Global nodal displacement field defining the prestress.
 * @return Mapped geometric stiffness in global translational DOF ordering.
 */
MapMatrix T3::stiffness_geom(Precision* buffer, const Field& displacement) {
    const Precision L0 = length_reference();
    logging::error(L0 > Precision(0),
        "T3: zero reference length in stiffness_geom for element ", this->elem_id);

    // Evaluate axial prestress from the supplied displacement field.
    const Vec3 u0 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
    const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[1]));

    const AxialStrainLinearized axial_strain((u1 - u0).dot(direction_reference()) / L0);
    AxialStressCauchy axial_stress;
    Precision         material_tangent = Precision(0);

    auto elasticity = get_elasticity();
    logging::error(elasticity->supports_axial_linearized(),
        "T3: material does not support linearized axial evaluation for element ", this->elem_id);

    // Read committed history without creating a persistent trial state.
    const Index      state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
    elasticity->evaluate(axial_strain, old_state, nullptr, axial_stress, material_tangent);

    // The geometric block acts isotropically on the translational difference
    // because the Total-Lagrangian truss stress term multiplies the identity.
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
 * Green-Lagrange strain is computed from the current stretch. The material reads
 * the committed state and writes the physical trial state. The resulting PK2
 * stress produces the internal force
 *
 *     f_int = A0 lambda S n
 *
 * and, when tangent storage is requested, the material and geometric blocks
 *
 *     K_mat = A0 C lambda^2/L0 (n tensor n),
 *     K_geo = A0 S/L0 I.
 *
 * Internal force is always assembled. A null matrix buffer therefore performs a
 * complete residual/material-state evaluation without constructing the tangent,
 * as required by residual-only nonlinear solver evaluations.
 *
 * @param buffer Optional six-by-six tangent storage; null requests force only.
 * @param nodal_forces Global nodal internal-force field to increment.
 * @param displacement Trial displacement field; current positions carry the
 *                     configuration used by the finite-strain truss kinematics.
 * @return Mapped tangent matrix, or an empty map when `buffer == nullptr`.
 */
MapMatrix T3::stiffness_tangent(Precision*   buffer,
                                NodeData&    nodal_forces,
                                const Field& displacement) {
    (void) displacement;

    logging::error(nodal_forces.components >= 3,
        "T3: nonlinear internal force requires at least three nodal components");

    // Evaluate stretch, Green-Lagrange strain and one physical constitutive
    // update from committed history into the persistent trial state.
    const Precision A0     = get_section()->area_;
    const Precision L0     = length_reference();
    const Precision lambda = stretch();

    logging::error(L0 > Precision(0),
        "T3: zero reference length in nonlinear tangent for element ", this->elem_id);

    const AxialStrainGreenLagrange strain = AxialStrainGreenLagrange::from_stretch(lambda);
    AxialStressPK2 stress;
    Precision      material_tangent = Precision(0);

    auto elasticity = get_elasticity();
    logging::error(elasticity->supports_axial_green_lagrange(),
        "T3: material does not support Green-Lagrange axial evaluation for element ", this->elem_id);

    const Index      state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
    Precision*       new_state = &(*this->_model_data->material_state_new)(state_row, 0);
    elasticity->evaluate(strain, old_state, new_state, stress, material_tangent);

    // Scatter the internal force directly from the constitutive stress. Stress
    // remains an element-local intermediate quantity throughout the evaluation.
    const Vec3  n      = direction_current();
    const Vec3  force  = A0 * lambda * stress.value() * n;
    const Index node_0 = static_cast<Index>(node_ids[0]);
    const Index node_1 = static_cast<Index>(node_ids[1]);

    for (Dim d = 0; d < 3; ++d) {
        nodal_forces(node_0, d) -= force(d);
        nodal_forces(node_1, d) += force(d);
    }

    // Residual-only evaluations stop after the physical material update and
    // internal-force assembly. A zero-sized map preserves the common return type
    // without requiring matrix storage from the caller.
    if (buffer == nullptr) {
        return MapMatrix(nullptr, 0, 0);
    }

    // Assemble material and stress-dependent tangent blocks from exactly the
    // same constitutive result used for the internal force.
    const Mat3 material_block =
        (A0 * material_tangent * lambda * lambda / L0) * (n * n.transpose());
    const Mat3 geometric_block = (A0 * stress.value() / L0) * Mat3::Identity();
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

Precision T3::integrate_scalar_field(bool scale_by_density, const ScalarField& field) {
    const Precision L = length_current();
    const Precision A = get_section()->area_;

    if (L <= Precision(0) || A <= Precision(0)) {
        return Precision(0);
    }

    return field(midpoint(*this)) * density_scale(*this, scale_by_density) * A * L;
}

Vec3 T3::integrate_vector_field(bool scale_by_density, const VecField& field) {
    const Precision L = length_current();
    const Precision A = get_section()->area_;

    if (L <= Precision(0) || A <= Precision(0)) {
        return Vec3::Zero();
    }

    return field(midpoint(*this)) * density_scale(*this, scale_by_density) * A * L;
}

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

Mat3 T3::integrate_tensor_field(bool scale_by_density, const TenField& field) {
    const Precision L = length_current();
    const Precision A = get_section()->area_;

    if (L <= Precision(0) || A <= Precision(0)) {
        return Mat3::Zero();
    }

    return field(midpoint(*this)) * density_scale(*this, scale_by_density) * A * L;
}

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
 * stress and pushes it forward as `sigma = lambda S` for physical output. The
 * element has one constant axial state, so the same values are copied to every
 * requested natural output coordinate.
 *
 * Result recovery reads committed material history but supplies no target state,
 * because recovery does not advance the nonlinear trial state.
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

    // Evaluate finite-strain PK2 stress and convert it to physical Cauchy stress.
    if (use_green_lagrange_nl) {
        const Precision lambda = stretch();
        const AxialStrainGreenLagrange axial_strain = AxialStrainGreenLagrange::from_stretch(lambda);
        AxialStressPK2 axial_stress;
        Precision      tangent = Precision(0);

        auto elasticity = get_elasticity();
        logging::error(elasticity->supports_axial_green_lagrange(),
            "T3: material does not support Green-Lagrange axial evaluation for element ", this->elem_id);
        elasticity->evaluate(axial_strain, old_state, nullptr, axial_stress, tangent);

        const AxialStressCauchy cauchy(lambda * axial_stress.value());
        strain_value = axial_strain.value();
        stress_value = cauchy.value();
    } else {
        // Evaluate infinitesimal axial strain and Cauchy stress on the reference axis.
        const Precision L0 = length_reference();
        logging::error(L0 > Precision(0),
            "T3: zero reference length in compute_stress_strain for element ", this->elem_id);

        const Vec3 u0 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
        const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[1]));

        const AxialStrainLinearized axial_strain((u1 - u0).dot(direction_reference()) / L0);
        AxialStressCauchy axial_stress;
        Precision         tangent = Precision(0);

        auto elasticity = get_elasticity();
        logging::error(elasticity->supports_axial_linearized(),
            "T3: material does not support linearized axial evaluation for element ", this->elem_id);
        elasticity->evaluate(axial_strain, old_state, nullptr, axial_stress, tangent);

        strain_value = axial_strain.value();
        stress_value = axial_stress.value();
    }

    // Replicate the constant axial state to all requested truss output points.
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
 * Relative nodal displacement is projected onto the reference truss direction
 * to obtain infinitesimal axial strain. The material is evaluated from committed
 * history without storing a target state, and Cauchy stress is multiplied by the
 * reference area. The constant axial force is copied into component zero at both
 * element nodes; all remaining section-force components are cleared.
 *
 * @param section_forces Element-nodal section-force output field.
 * @param displacement Global nodal displacement field.
 * @param offset First element-nodal output row belonging to this truss.
 * @return Always `true` after both nodal result rows were written.
 */
bool T3::compute_beam_section_forces(Field&       section_forces,
                                     const Field& displacement,
                                     int          offset) {
    const Vec3 u0 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
    const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[1]));

    const Precision L0 = length_reference();
    logging::error(L0 > Precision(0),
        "T3: zero reference length in compute_beam_section_forces for element ", this->elem_id);

    const AxialStrainLinearized axial_strain((u1 - u0).dot(direction_reference()) / L0);
    AxialStressCauchy axial_stress;
    Precision         tangent = Precision(0);

    auto elasticity = get_elasticity();
    logging::error(elasticity->supports_axial_linearized(),
        "T3: material does not support linearized axial evaluation for element ", this->elem_id);

    const Index      state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
    elasticity->evaluate(axial_strain, old_state, nullptr, axial_stress, tangent);

    const Precision axial_force = get_section()->area_ * axial_stress.value();

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
