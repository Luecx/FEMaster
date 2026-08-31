/**
 * @file truss.cpp
 * @brief Implements the three-dimensional two-node truss element.
 *
 * The nonlinear Total-Lagrangian path performs exactly one source-to-target
 * constitutive integration per material point and solver evaluation. Read-only
 * stiffness and result-recovery paths evaluate the immutable source history and
 * can therefore never advance plastic state.
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
    if (!scale_by_density) return Precision(1);

    auto mat = elem.get_material();
    logging::error(mat && mat->has_density(),
        "T3: material density is required when scale_by_density=true for element ", elem.elem_id);
    return mat->get_density();
}

const Precision* source_state(T3& elem, material::Elasticity* elasticity) {
    if (elasticity->state_size() == 0) return nullptr;

    logging::error(elem._model_data->material_state != nullptr,
        "T3: stateful material requires source material state for element ", elem.elem_id);
    return &(*elem._model_data->material_state)(elem.mp_index(0), 0);
}

Precision* target_state(T3& elem, material::Elasticity* elasticity) {
    if (elasticity->state_size() == 0) return nullptr;

    logging::error(elem._model_data->material_state_target != nullptr,
        "T3: stateful material requires target material state for element ", elem.elem_id);
    return &(*elem._model_data->material_state_target)(elem.mp_index(0), 0);
}

} // namespace

T3::T3(ID elem_id, std::array<ID, N> node_ids_in)
    : StructuralElement(elem_id), node_ids(node_ids_in) {}

ElDofs T3::dofs() const { return ElDofs{true, true, true, false, false, false}; }
Dim T3::dimensions() const { return 3; }
Dim T3::n_nodes() const { return N; }
Dim T3::num_ip() const { return 1; }
const ID* T3::nodes() const { return node_ids.data(); }

SurfacePtr T3::surface(ID surface_id) {
    (void) surface_id;
    return nullptr;
}

std::string T3::type_name() const { return "T3"; }

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
    logging::error(this->_model_data != nullptr && this->_model_data->positions_reference != nullptr,
        "T3: reference positions are not available for element ", this->elem_id);
    return this->_model_data->positions_reference->row_vec3(static_cast<Index>(node_ids[local_node]));
}

Vec3 T3::node_position_current(Index local_node) const {
    logging::error(local_node < N,
        "T3: local node index out of range in element ", this->elem_id);
    logging::error(this->_model_data != nullptr && this->_model_data->positions != nullptr,
        "T3: current positions are not available for element ", this->elem_id);
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
    logging::error(L0 > Precision(0), "T3: zero reference length for element ", elem_id);
    logging::error(l > Precision(0), "T3: zero current length for element ", elem_id);
    return l / L0;
}

Vec3 T3::direction_reference() const {
    const Precision L0 = length_reference();
    logging::error(L0 > Precision(0), "T3: zero reference length in element ", elem_id);
    return (node_position_reference(1) - node_position_reference(0)) / L0;
}

Vec3 T3::direction_current() const {
    const Precision l = length_current();
    logging::error(l > Precision(0), "T3: zero current length in element ", elem_id);
    return (node_position_current(1) - node_position_current(0)) / l;
}

Precision T3::length() { return length_current(); }
Vec3 T3::direction() { return direction_current(); }
Precision T3::volume() { return get_section()->area_ * length_current(); }

/**
 * Assembles state-neutral Total-Lagrangian material stiffness.
 */
MapMatrix T3::stiffness(Precision* buffer) {
    const Precision A0     = get_section()->area_;
    const Precision L0     = length_reference();
    const Precision lambda = stretch();
    logging::error(L0 > Precision(0), "T3: zero reference length in stiffness for element ", elem_id);

    const AxialStrainGreenLagrange strain = AxialStrainGreenLagrange::from_stretch(lambda);
    AxialStressPK2 stress;
    Precision material_tangent = Precision(0);

    auto elasticity = get_elasticity();
    logging::error(elasticity->supports_axial_green_lagrange(),
        "T3: material does not support Green-Lagrange axial evaluation for element ", elem_id);

    elasticity->evaluate(strain, source_state(*this, elasticity), stress, material_tangent);

    const Vec3 n = direction_current();
    const Mat3 k = (A0 * material_tangent * lambda * lambda / L0) * (n * n.transpose());

    StaticMatrix<N * 3, N * 3> K = StaticMatrix<N * 3, N * 3>::Zero();
    K.block(0, 0, 3, 3) =  k;
    K.block(0, 3, 3, 3) = -k;
    K.block(3, 0, 3, 3) = -k;
    K.block(3, 3, 3, 3) =  k;

    MapMatrix result(buffer, N * 3, N * 3);
    result = K;
    return result;
}

MapMatrix T3::stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) {
    logging::error(ip_stress.components >= 1,
        "T3: geometric stiffness requires nonlinear IP stress component 0");

    const Precision L0 = length_reference();
    logging::error(L0 > Precision(0), "T3: zero reference length in stiffness_geom for element ", elem_id);

    const Precision second_piola = ip_stress(static_cast<Index>(ip_start_idx), 0);
    const Mat3 k = (get_section()->area_ * second_piola / L0) * Mat3::Identity();

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
 * Assembles tangent and internal force from one explicit constitutive integration.
 * Source history is immutable; the complete updated state is written to the
 * model's target material field.
 */
MapMatrix T3::stiffness_tangent(Precision*   buffer,
                                Field&       ip_stress_state,
                                NodeData&    nodal_forces,
                                const Field& displacement) {
    (void) displacement;

    logging::error(ip_stress_state.components >= 1,
        "T3: nonlinear stress state requires at least one component");
    logging::error(nodal_forces.components >= 3,
        "T3: nonlinear internal force requires at least three nodal components");

    const Precision A0     = get_section()->area_;
    const Precision L0     = length_reference();
    const Precision lambda = stretch();
    logging::error(L0 > Precision(0),
        "T3: zero reference length in nonlinear tangent for element ", elem_id);

    const AxialStrainGreenLagrange strain = AxialStrainGreenLagrange::from_stretch(lambda);
    AxialStressPK2 stress;
    Precision material_tangent = Precision(0);

    auto elasticity = get_elasticity();
    logging::error(elasticity->supports_axial_green_lagrange(),
        "T3: material does not support Green-Lagrange axial integration for element ", elem_id);

    elasticity->integrate(
        strain,
        source_state(*this, elasticity),
        target_state(*this, elasticity),
        stress,
        material_tangent
    );

    ip_stress_state(this->ip_index(0), 0) = stress.value();

    const Vec3 n = direction_current();
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

    const Vec3 force = A0 * lambda * stress.value() * n;
    const Index node_0 = static_cast<Index>(node_ids[0]);
    const Index node_1 = static_cast<Index>(node_ids[1]);

    for (Dim d = 0; d < 3; ++d) {
        nodal_forces(node_0, d) -= force(d);
        nodal_forces(node_1, d) += force(d);
    }

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
    if (L <= Precision(0) || A <= Precision(0)) return Precision(0);
    return field(midpoint(*this)) * density_scale(*this, scale_by_density) * A * L;
}

Vec3 T3::integrate_vector_field(bool scale_by_density, const VecField& field) {
    const Precision L = length_current();
    const Precision A = get_section()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return Vec3::Zero();
    return field(midpoint(*this)) * density_scale(*this, scale_by_density) * A * L;
}

void T3::integrate_vector_field(Field& node_loads, bool scale_by_density, const VecField& field) {
    const Precision L = length_current();
    const Precision A = get_section()->area_;
    if (L <= Precision(0) || A <= Precision(0)) return;

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
    if (L <= Precision(0) || A <= Precision(0)) return Mat3::Zero();
    return field(midpoint(*this)) * density_scale(*this, scale_by_density) * A * L;
}

void T3::apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) {
    (void) node_loads;
    (void) node_temp;
    (void) ref_temp;
}

/**
 * Recovers strain and physical stress without advancing constitutive history.
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
    auto elasticity = get_elasticity();
    const Precision* state = source_state(*this, elasticity);

    if (use_green_lagrange_nl) {
        const Precision lambda = stretch();
        const AxialStrainGreenLagrange axial_strain =
            AxialStrainGreenLagrange::from_stretch(lambda);
        AxialStressPK2 axial_stress;
        Precision tangent = Precision(0);

        logging::error(elasticity->supports_axial_green_lagrange(),
            "T3: material does not support Green-Lagrange axial evaluation for element ", elem_id);
        elasticity->evaluate(axial_strain, state, axial_stress, tangent);

        strain_value = axial_strain.value();
        stress_value = lambda * axial_stress.value();
    } else {
        const Precision L0 = length_reference();
        logging::error(L0 > Precision(0),
            "T3: zero reference length in compute_stress_strain for element ", elem_id);

        const Vec3 u0 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
        const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[1]));
        const AxialStrainLinearized axial_strain((u1 - u0).dot(direction_reference()) / L0);
        AxialStressCauchy axial_stress;
        Precision tangent = Precision(0);

        logging::error(elasticity->supports_axial_linearized(),
            "T3: material does not support linearized axial evaluation for element ", elem_id);
        elasticity->evaluate(axial_strain, state, axial_stress, tangent);

        strain_value = axial_strain.value();
        stress_value = axial_stress.value();
    }

    for (Index i = 0; i < static_cast<Index>(rst.rows()); ++i) {
        const Index row = static_cast<Index>(offset) + i;

        if (strain) {
            for (Index j = 0; j < strain->components; ++j) (*strain)(row, j) = Precision(0);
            (*strain)(row, 0) = strain_value;
        }

        if (stress) {
            for (Index j = 0; j < stress->components; ++j) (*stress)(row, j) = Precision(0);
            (*stress)(row, 0) = stress_value;
        }
    }
}

void T3::compute_stress_state(Field&       stress_state,
                              const Field& displacement,
                              int          offset,
                              bool         use_green_lagrange_nl) {
    if (use_green_lagrange_nl) {
        const Precision lambda = stretch();
        const AxialStrainGreenLagrange strain = AxialStrainGreenLagrange::from_stretch(lambda);
        AxialStressPK2 stress;
        Precision tangent = Precision(0);

        auto elasticity = get_elasticity();
        logging::error(elasticity->supports_axial_green_lagrange(),
            "T3: material does not support Green-Lagrange axial evaluation for element ", elem_id);
        elasticity->evaluate(strain, source_state(*this, elasticity), stress, tangent);

        const RowMatrix rst = stress_strain_ip_rst();
        for (Index i = 0; i < static_cast<Index>(rst.rows()); ++i) {
            const Index row = static_cast<Index>(offset) + i;
            for (Index component = 0; component < stress_state.components; ++component) {
                stress_state(row, component) = Precision(0);
            }
            stress_state(row, 0) = stress.value();
        }
        return;
    }

    compute_stress_strain(nullptr, &stress_state, displacement,
                          stress_strain_ip_rst(), offset, false);
}

void T3::compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) {
    logging::error(node_forces.domain == FieldDomain::NODE,
        "T3: nonlinear internal force output must use NODE domain");
    logging::error(node_forces.components >= 3,
        "T3: nonlinear internal force output requires at least 3 components");
    logging::error(ip_stress.components >= 1,
        "T3: nonlinear internal force requires IP stress component 0");

    const Precision lambda        = stretch();
    const Precision second_piola = ip_stress(ip_index(0), 0);
    const Vec3 force = get_section()->area_ * lambda * second_piola * direction_current();

    const Index node_0 = static_cast<Index>(node_ids[0]);
    const Index node_1 = static_cast<Index>(node_ids[1]);
    for (Dim d = 0; d < 3; ++d) {
        node_forces(node_0, d) -= force(d);
        node_forces(node_1, d) += force(d);
    }
}

void T3::compute_compliance(Field& displacement, Field& result) {
    Precision buffer[N * 3 * N * 3] {};
    MapMatrix K = stiffness(buffer);

    StaticVector<N * 3> u;
    for (Index i = 0; i < N; ++i) {
        const Vec3 ui = displacement.row_vec3(static_cast<Index>(node_ids[i]));
        for (Index d = 0; d < 3; ++d) u(i * 3 + d) = ui(d);
    }

    result(static_cast<Index>(this->elem_id), 0) = u.dot(K * u);
}

bool T3::compute_beam_section_forces(Field&       section_forces,
                                     const Field& displacement,
                                     int          offset) {
    const Vec3 u0 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
    const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[1]));

    const Precision L0 = length_reference();
    logging::error(L0 > Precision(0),
        "T3: zero reference length in compute_beam_section_forces for element ", elem_id);

    const AxialStrainLinearized strain((u1 - u0).dot(direction_reference()) / L0);
    AxialStressCauchy stress;
    Precision tangent = Precision(0);

    auto elasticity = get_elasticity();
    logging::error(elasticity->supports_axial_linearized(),
        "T3: material does not support linearized axial evaluation for element ", elem_id);
    elasticity->evaluate(strain, source_state(*this, elasticity), stress, tangent);

    const Precision axial_force = get_section()->area_ * stress.value();
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
