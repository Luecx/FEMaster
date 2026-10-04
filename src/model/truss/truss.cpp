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

#include "../../material/isotropic_j2_elasticity.h"

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
 * Evaluates the truss response about a selected linearization state u0.
 *
 * A null linearization denotes u0 = 0. The complete tangent is evaluated at u0
 * and contains both the material and geometric contributions generated by the
 * stress already present at that state. Internal force away from u0 is returned
 * as the affine approximation
 *
 *     f(u) = f(u0) + K_T(u0) (u - u0).
 *
 * The separately requested geometric stiffness is generated only by the
 * linearized stress increment from u0 to the requested state u. Existing stress
 * at u0 therefore remains part of K_T(u0) and is not repeated in the separate
 * geometric operator.
 *
 * Trial constitutive history may only be written for an exact evaluation at u0.
 */
MapMatrix T3::evaluate(
    Precision*   tangent_buffer,
    Precision*   geometric_tangent_buffer,
    NodeData*    internal_force,
    const Field* displacement,
    const Field* linearization,
    const Field* thermal_free_strain,
    bool         update_state
) {
    (void) thermal_free_strain;


    // -----------------------------------------------------------------------------
    // requested outputs
    // -----------------------------------------------------------------------------

    // Requested element outputs:
    // - complete tangent stiffness at the linearization state u0
    // - geometric stiffness generated by the stress increment from u0 to u
    // - internal force at the requested displacement state u
    const bool with_tangent   = tangent_buffer           != nullptr;
    const bool with_geometric = geometric_tangent_buffer != nullptr;
    const bool with_force     = internal_force           != nullptr;


    // -----------------------------------------------------------------------------
    // reference geometry and evaluation state
    // -----------------------------------------------------------------------------

    // Reference section area and reference element length used by both the
    // linearized and Total-Lagrangian truss formulations.
    const Precision A0 = get_section()->area_;
    const Precision L0 = length_reference();

    logging::error(with_tangent || with_geometric || with_force,
        "T3: evaluation requires at least one requested output");
    logging::error(!with_force || displacement != nullptr,
        "T3: internal force evaluation requires displacement");
    logging::error(!update_state || (linearization != nullptr && displacement == linearization),
        "T3: material state requires an exact evaluation at the linearization state");
    logging::error(L0 > Precision(0),
        "T3: zero reference length for element ", this->elem_id);

    // Gather the displacement at the linearization state u0 and the displacement
    // increment Delta u = u - u0 for both element nodes.
    StaticMatrix<N, 3> disp_linear_base = StaticMatrix<N, 3>::Zero();
    StaticMatrix<N, 3> disp_delta       = StaticMatrix<N, 3>::Zero();

    if (linearization) {
        disp_linear_base.row(0) = linearization->row_vec3(static_cast<Index>(node_ids[0])).transpose();
        disp_linear_base.row(1) = linearization->row_vec3(static_cast<Index>(node_ids[1])).transpose();
    }

    if (displacement) {
        disp_delta.row(0) = displacement->row_vec3(static_cast<Index>(node_ids[0])).transpose() - disp_linear_base.row(0);
        disp_delta.row(1) = displacement->row_vec3(static_cast<Index>(node_ids[1])).transpose() - disp_linear_base.row(1);
    }

    const Vec3 delta_axis = (disp_delta.row(1) - disp_delta.row(0)).transpose();


    // -----------------------------------------------------------------------------
    // required internal quantities
    // -----------------------------------------------------------------------------

    // An affine force evaluation is required whenever the requested state u
    // differs from the linearization state u0:
    //
    //     f(u) = f(u0) + K_T(u0) (u - u0).
    const bool affine_force = with_force && displacement != linearization;

    // The complete tangent is needed either as an explicit output or internally
    // for the affine force correction.
    const bool need_tangent = with_tangent || affine_force;

    // The constitutive tangent is additionally required to obtain the stress
    // increment used by the separately requested geometric stiffness.
    const bool need_material = need_tangent || with_geometric;


    // -----------------------------------------------------------------------------
    // material state
    // -----------------------------------------------------------------------------

    auto elasticity = get_elasticity();

    // Locate the single constitutive state belonging to this truss. old_state
    // contains the committed history; new_state is writable only for an exact
    // physical evaluation at u0.
    const Index      state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
    Precision*       new_state = update_state ? &(*this->_model_data->material_state_new)(state_row, 0) : nullptr;


    // -----------------------------------------------------------------------------
    // kinematics and constitutive response at u0
    // -----------------------------------------------------------------------------

    // The quantities below describe the linearization state u0 and are used by
    // all subsequent force and stiffness calculations.
    //
    // For u0 = 0 the truss uses infinitesimal axial strain on the reference
    // axis. For a finite u0 it uses Total-Lagrangian Green-Lagrange strain and
    // the corresponding PK2 stress.
    Vec3      direction_base  = direction_reference();
    Precision stretch_base    = Precision(1);
    Precision stress_base     = Precision(0);
    Precision material_tangent = Precision(0);

    if (linearization == nullptr) {
        logging::error(elasticity->supports_axial_linearized(),
            "T3: material does not support linearized axial evaluation for element ",
            this->elem_id);

        AxialStressCauchy stress;
        elasticity->evaluate(
            AxialStrainLinearized(Precision(0)),
            old_state,
            new_state,
            stress,
            need_material ? &material_tangent : nullptr
        );

        stress_base = stress.value();

    } else {
        logging::error(elasticity->supports_axial_green_lagrange(),
            "T3: material does not support Green-Lagrange axial evaluation for element ",
            this->elem_id);

        // Construct the element axis in the linearization configuration
        //
        //     x0 = X + u0.
        const Vec3 X0 = node_position_reference(0);
        const Vec3 X1 = node_position_reference(1);
        const Vec3 x0 = X0 + disp_linear_base.row(0).transpose();
        const Vec3 x1 = X1 + disp_linear_base.row(1).transpose();
        const Vec3 axis_base = x1 - x0;

        const Precision length_base = axis_base.norm();
        logging::error(length_base > Precision(0),
            "T3: zero length at the linearization state for element ", this->elem_id);

        stretch_base   = length_base / L0;
        direction_base = axis_base / length_base;

        // Green-Lagrange axial strain at the linearization state:
        //
        //     E0 = 1/2 (lambda0^2 - 1).
        const AxialStrainGreenLagrange strain =
            AxialStrainGreenLagrange::from_stretch(stretch_base);

        AxialStressPK2 stress;
        elasticity->evaluate(
            strain,
            old_state,
            new_state,
            stress,
            need_material ? &material_tangent : nullptr
        );

        stress_base = stress.value();
    }


    // -----------------------------------------------------------------------------
    // complete tangent stiffness at u0
    // -----------------------------------------------------------------------------

    StaticMatrix<N * 3, N * 3> tangent = StaticMatrix<N * 3, N * 3>::Zero();

    if (need_tangent) {
        // The consistent Total-Lagrangian tangent consists of
        //
        //     K_T(u0) = K_M(u0) + K_G(S0).
        //
        // At the reference state stretch_base = 1, so the same expression
        // reduces directly to the ordinary linear truss stiffness.
        const Mat3 material_block  = (A0 * material_tangent * stretch_base * stretch_base / L0)
                                   * (direction_base * direction_base.transpose());
        const Mat3 geometric_block = (A0 * stress_base / L0) * Mat3::Identity();
        const Mat3 tangent_block   = material_block + geometric_block;

        tangent.block(0, 0, 3, 3) =  tangent_block;
        tangent.block(0, 3, 3, 3) = -tangent_block;
        tangent.block(3, 0, 3, 3) = -tangent_block;
        tangent.block(3, 3, 3, 3) =  tangent_block;
    }


    // -----------------------------------------------------------------------------
    // geometric stiffness generated by the perturbation from u0 to u
    // -----------------------------------------------------------------------------

    StaticMatrix<N * 3, N * 3> geometric = StaticMatrix<N * 3, N * 3>::Zero();

    if (with_geometric) {
        // Linearize the Green-Lagrange axial strain about u0:
        //
        //     Delta E = lambda0 n0 . Delta(x1 - x0) / L0.
        //
        // For u0 = 0, lambda0 = 1 and n0 is the reference direction, so this
        // reduces to the ordinary infinitesimal axial strain increment.
        const Precision strain_increment = stretch_base * direction_base.dot(delta_axis) / L0;

        // Convert the strain increment into the corresponding stress increment
        // using the constitutive tangent evaluated at u0:
        //
        //     Delta S = C0 Delta E.
        //
        // Only this stress increment enters the separately requested geometric
        // stiffness. The base stress S0 is already contained in K_T(u0).
        const Precision stress_increment = material_tangent * strain_increment;
        const Mat3 geometric_block       = (A0 * stress_increment / L0) * Mat3::Identity();

        geometric.block(0, 0, 3, 3) =  geometric_block;
        geometric.block(0, 3, 3, 3) = -geometric_block;
        geometric.block(3, 0, 3, 3) = -geometric_block;
        geometric.block(3, 3, 3, 3) =  geometric_block;
    }


    // -----------------------------------------------------------------------------
    // internal force
    // -----------------------------------------------------------------------------

    if (with_force) {
        StaticVector<N * 3> force = StaticVector<N * 3>::Zero();

        // Internal force at the linearization state. In the Total-Lagrangian
        // formulation the reference-area PK2 stress produces the axial force
        //
        //     N0 = A0 lambda0 S0.
        //
        // For u0 = 0 this reduces to A0 S0.
        const Vec3 axial_force = A0 * stretch_base * stress_base * direction_base;

        force.template segment<3>(0) = -axial_force;
        force.template segment<3>(3) =  axial_force;

        // If u differs from u0, continue the force affinely with the complete
        // tangent evaluated at u0.
        if (affine_force) {
            StaticVector<N * 3> delta = StaticVector<N * 3>::Zero();
            delta.template segment<3>(0) = disp_delta.row(0).transpose();
            delta.template segment<3>(3) = disp_delta.row(1).transpose();

            force.noalias() += tangent * delta;
        }

        const Index node0 = static_cast<Index>(node_ids[0]);
        const Index node1 = static_cast<Index>(node_ids[1]);

        for (Dim d = 0; d < 3; ++d) {
            (*internal_force)(node0, d) += force(d);
            (*internal_force)(node1, d) += force(3 + d);
        }
    }


    // -----------------------------------------------------------------------------
    // copy requested operators into caller-owned storage
    // -----------------------------------------------------------------------------

    if (with_tangent) {
        MapMatrix mapped(tangent_buffer, N * 3, N * 3);
        mapped = tangent;
    }

    if (with_geometric) {
        MapMatrix mapped(geometric_tangent_buffer, N * 3, N * 3);
        mapped = geometric;
    }

    return with_tangent   ? MapMatrix(tangent_buffer,           N * 3, N * 3)
         : with_geometric ? MapMatrix(geometric_tangent_buffer, N * 3, N * 3)
                          : MapMatrix(nullptr, 0, 0);
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
 * @param linearization Null for reference recovery, displacement for exact recovery.
 * @param thermal_free_strain Unused; axial thermal recovery remains unchanged.
 */
void T3::compute_stress_strain(
    Field*           strain,
    Field*           stress,
    const Field&     displacement,
    const RowMatrix& rst,
    int              offset,
    const Field*     linearization,
    const Field*     thermal_free_strain
) {
    // This formulation supports reference and exact recovery states
    logging::error(linearization == nullptr || linearization == &displacement,
        "T3: intermediate recovery expansion points are not supported");
    const bool use_green_lagrange_nl = linearization != nullptr;
    (void) thermal_free_strain;

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
 * Recovers accumulated equivalent plastic strain from the committed truss state.
 *
 * The truss owns one constitutive material point, so its scalar PEEQ value is
 * constant over the element and is copied to both element-nodal rows. Recovery
 * reads accepted J2 history directly and does not reevaluate the material law.
 *
 * A J2 material whose nonlinear state storage has not been initialized yet
 * contributes zero. Non-J2 materials return false and are excluded from the
 * model-wide nodal average.
 *
 * @param peeq Scalar ELEMENT_NODAL output field.
 * @param offset First element-nodal row belonging to the truss.
 * @return True when the truss material uses J2 plasticity.
 */
bool T3::compute_peeq(Field& peeq, int offset) {
    logging::error(peeq.domain == FieldDomain::ELEMENT_NODAL && peeq.components == 1,
        "T3: PEEQ recovery requires scalar ELEMENT_NODAL output");

    auto mat = get_material();
    if (!mat || !mat->has_elasticity()) return false;

    const auto* j2 = mat->elasticity()->as<material::IsotropicJ2Elasticity>();
    if (!j2) return false;

    Precision value = Precision(0);
    const auto& state = this->_model_data->material_state_old;
    if (state && state->components >= j2->state_size()) {
        const Precision* old_state = &(*state)(this->mp_index(0), 0);
        value = j2->equivalent_plastic_strain(old_state);
    }

    peeq(static_cast<Index>(offset) + 0, 0) = value;
    peeq(static_cast<Index>(offset) + 1, 0) = value;
    return true;
}

/**
 * Computes the element compliance contribution `u^T K u`.
 *
 * @param displacement Global nodal displacement field.
 * @param result Element result field receiving the scalar compliance.
 */
void T3::compute_compliance(Field& displacement, Field& result) {
    Precision buffer[N * 3 * N * 3] {};
    MapMatrix K = evaluate(buffer, nullptr, nullptr, nullptr, nullptr, nullptr, false);

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