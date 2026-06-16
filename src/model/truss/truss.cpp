#include "truss.h"

#include <cmath>

namespace fem {
namespace model {
namespace {

Vec3 midpoint(T3& elem) {
    return (elem.node_position(0) + elem.node_position(1)) * Precision(0.5);
}

Precision density_scale(T3& elem, bool scale_by_density) {
    if (!scale_by_density) {
        return Precision(1);
    }

    auto mat = elem.get_material();
    logging::error(mat && mat->has_density(),
                   "T3: material density is required when scale_by_density=true for element ",
                   elem.elem_id);
    return mat->get_density();
}

} // namespace

T3::T3(ID elem_id, std::array<ID, N> node_ids_in)
    : StructuralElement(elem_id)
    , node_ids(node_ids_in) {}

ElDofs T3::dofs() {
    return ElDofs{true, true, true, false, false, false};
}

Dim T3::dimensions() {
    return 3;
}

Dim T3::n_nodes() {
    return N;
}

Dim T3::n_integration_points() {
    return 1;
}

ID* T3::nodes() {
    return node_ids.data();
}

SurfacePtr T3::surface(ID surface_id) {
    (void)surface_id;
    return nullptr;
}

std::string T3::type_name() const {
    return "T3";
}

TrussSection* T3::get_section() {
    logging::error(this->_section != nullptr,
                   "T3: missing section for element ",
                   this->elem_id);
    auto* section = this->_section->template as<TrussSection>();
    logging::error(section != nullptr,
                   "T3: section is not a truss section for element ",
                   this->elem_id);
    return section;
}

material::MaterialPtr T3::get_material() {
    TrussSection* section = get_section();
    logging::error(section->material_ != nullptr,
                   "T3: no material set for element ",
                   this->elem_id);
    return section->material_;
}

material::IsotropicElasticity* T3::get_elasticity() {
    auto mat = get_material();
    logging::error(mat->has_elasticity(),
                   "T3: material has no elasticity for element ",
                   this->elem_id);
    auto* elasticity = mat->elasticity()->template as<material::IsotropicElasticity>();
    logging::error(elasticity != nullptr,
                   "T3: material elasticity is not isotropic for element ",
                   this->elem_id);
    return elasticity;
}

Precision T3::length() {
    return (this->node_position(1) - this->node_position(0)).norm();
}

Vec3 T3::direction() {
    const Precision L = length();
    logging::error(L > Precision(0),
                   "T3: zero length for element ",
                   this->elem_id);
    return (this->node_position(1) - this->node_position(0)) / L;
}

Precision T3::volume() {
    return integrate_scalar_field(false, [](const Vec3&) { return Precision(1); });
}

MapMatrix T3::stiffness(Precision* buffer) {
    const Precision E = get_elasticity()->youngs;
    const Precision A = get_section()->area_;
    const Precision L = length();
    logging::error(L > Precision(0),
                   "T3: zero length in stiffness for element ",
                   this->elem_id);

    StaticMatrix<2, 2> k_local;
    k_local << Precision(1), Precision(-1),
               Precision(-1), Precision(1);

    const Vec3 t = direction();
    StaticMatrix<N * 3, 2> P = StaticMatrix<N * 3, 2>::Zero();
    P.block(0, 0, 3, 1) = t;
    P.block(3, 1, 3, 1) = t;

    MapMatrix result(buffer, N * 3, N * 3);
    result = P * (k_local * (E * A / L)) * P.transpose();
    return result;
}

MapMatrix T3::stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) {
    logging::error(ip_stress.components >= 1,
                   "T3: geometric stiffness requires nonlinear IP stress component 0");

    const Precision L = length();
    logging::error(L > Precision(0),
                   "T3: zero length in stiffness_geom for element ",
                   this->elem_id);

    const Vec3 n = direction();
    const Precision sigma = ip_stress(static_cast<Index>(ip_start_idx), 0);
    const Precision axial_force = get_section()->area_ * sigma;
    const Mat3 projector = Mat3::Identity() - n * n.transpose();
    const Mat3 block = (axial_force / L) * projector;

    StaticMatrix<N * 3, N * 3> Kg = StaticMatrix<N * 3, N * 3>::Zero();
    Kg.block(0, 0, 3, 3) = block;
    Kg.block(0, 3, 3, 3) = -block;
    Kg.block(3, 0, 3, 3) = -block;
    Kg.block(3, 3, 3, 3) = block;

    MapMatrix result(buffer, N * 3, N * 3);
    result = Kg;
    return result;
}

MapMatrix T3::mass(Precision* buffer) {
    StaticMatrix<N * 3, N * 3> M = StaticMatrix<N * 3, N * 3>::Zero();

    auto mat = get_material();
    if (mat->has_density()) {
        const Precision rho = mat->get_density();
        const Precision A = get_section()->area_;
        const Precision L = length();
        const Precision m = rho * A * L;

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

Precision T3::integrate_scalar_field(bool scale_by_density,
                                     const ScalarField& field) {
    const Precision L = length();
    const Precision A = get_section()->area_;
    if (L <= Precision(0) || A <= Precision(0)) {
        return Precision(0);
    }
    return field(midpoint(*this)) * (density_scale(*this, scale_by_density) * A * L);
}

Vec3 T3::integrate_vector_field(bool scale_by_density,
                                const VecField& field) {
    const Precision L = length();
    const Precision A = get_section()->area_;
    if (L <= Precision(0) || A <= Precision(0)) {
        return Vec3::Zero();
    }
    return field(midpoint(*this)) * (density_scale(*this, scale_by_density) * A * L);
}

void T3::integrate_vector_field(Field& node_loads,
                                bool scale_by_density,
                                const VecField& field) {
    const Precision L = length();
    const Precision A = get_section()->area_;
    if (L <= Precision(0) || A <= Precision(0)) {
        return;
    }

    const Vec3 F = field(midpoint(*this)) * (density_scale(*this, scale_by_density) * A * L);
    for (Index i = 0; i < N; ++i) {
        const Index n_id = static_cast<Index>(node_ids[i]);
        node_loads(n_id, 0) += F(0) * Precision(0.5);
        node_loads(n_id, 1) += F(1) * Precision(0.5);
        node_loads(n_id, 2) += F(2) * Precision(0.5);
    }
}

Mat3 T3::integrate_tensor_field(bool scale_by_density,
                                const TenField& field) {
    const Precision L = length();
    const Precision A = get_section()->area_;
    if (L <= Precision(0) || A <= Precision(0)) {
        return Mat3::Zero();
    }
    return field(midpoint(*this)) * (density_scale(*this, scale_by_density) * A * L);
}

void T3::apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) {
    (void)node_loads;
    (void)node_temp;
    (void)ref_temp;
}

void T3::compute_stress_strain(Field* strain,
                               Field* stress,
                               const Field& displacement,
                               const RowMatrix& rst,
                               int offset,
                               bool use_green_lagrange_nl) {
    logging::error(strain != nullptr || stress != nullptr,
                   "T3: compute_stress_strain requires at least one output field");
    logging::error(rst.cols() >= 1,
                   "T3: stress/strain coordinates require at least 1 column");

    Precision strain_value = Precision(0);
    Precision stress_value = Precision(0);

    if (use_green_lagrange_nl) {
        const Vec3 X0 = this->node_position(0);
        const Vec3 X1 = this->node_position(1);
        const Vec3 u0 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
        const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[1]));
        const Vec3 x0 = X0 + u0;
        const Vec3 x1 = X1 + u1;

        const Precision L0 = (X1 - X0).norm();
        const Precision l = (x1 - x0).norm();
        logging::error(L0 > Precision(0),
                       "T3: zero reference length in compute_stress_strain for element ",
                       this->elem_id);
        logging::error(l > Precision(0),
                       "T3: zero current length in compute_stress_strain for element ",
                       this->elem_id);

        const Precision lambda = l / L0;
        strain_value = Precision(0.5) * (lambda * lambda - Precision(1));
        const Precision second_piola_stress = get_elasticity()->youngs * strain_value;
        stress_value = lambda * second_piola_stress;
    } else {
        const Precision L = length();
        logging::error(L > Precision(0),
                       "T3: zero length in compute_stress_strain for element ",
                       this->elem_id);
        const Vec3 u0 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
        const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[1]));
        strain_value = (u1 - u0).dot(direction()) / L;
        stress_value = get_elasticity()->youngs * strain_value;
    }

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

void T3::compute_stress_state(Field& stress_state,
                              const Field& displacement,
                              int offset,
                              bool use_green_lagrange_nl) {
    RowMatrix rst = stress_strain_ip_rst();
    compute_stress_strain(nullptr, &stress_state, displacement, rst, offset, use_green_lagrange_nl);
}

void T3::compute_internal_force_nonlinear(Field& node_forces,
                                          const Field& displacement,
                                          const Field& ip_stress,
                                          int ip_offset) {
    (void) displacement;

    logging::error(node_forces.domain == FieldDomain::NODE,
                   "T3: nonlinear internal force output must use NODE domain");
    logging::error(node_forces.components >= 3,
                   "T3: nonlinear internal force output requires at least 3 components");
    logging::error(ip_stress.components >= 1,
                   "T3: nonlinear internal force requires IP stress component 0");

    const Precision L = length();
    logging::error(L > Precision(0),
                   "T3: zero length in compute_internal_force_nonlinear for element ",
                   this->elem_id);

    const Vec3 n = direction();
    const Precision sigma = ip_stress(static_cast<Index>(ip_offset), 0);
    const Vec3 force = get_section()->area_ * sigma * n;

    const Index n0 = static_cast<Index>(node_ids[0]);
    const Index n1 = static_cast<Index>(node_ids[1]);
    for (Dim d = 0; d < 3; ++d) {
        node_forces(n0, d) -= force(d);
        node_forces(n1, d) += force(d);
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

void T3::compute_compliance_angle_derivative(Field& displacement, Field& result) {
    (void)displacement;
    (void)result;
}

bool T3::compute_shear_flow(Field& shear_flow,
                            const Field& displacement,
                            int offset) {
    (void)shear_flow;
    (void)displacement;
    (void)offset;
    return false;
}

bool T3::compute_beam_section_forces(Field& section_forces,
                                     const Field& displacement,
                                     int offset) {
    const Vec3 u0 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
    const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[1]));

    const Precision L = length();
    logging::error(L > Precision(0),
                   "T3: zero length in compute_beam_section_forces for element ",
                   this->elem_id);

    const Precision axial_strain = (u1 - u0).dot(direction()) / L;
    const Precision axial_force = get_elasticity()->youngs * get_section()->area_ * axial_strain;

    for (Index i = 0; i < N; ++i) {
        for (Index d = 0; d < section_forces.components; ++d) {
            section_forces(static_cast<Index>(offset) + i, d) = Precision(0);
        }
        section_forces(static_cast<Index>(offset) + i, 0) = axial_force;
    }
    return true;
}

bool T3::compute_shell_section_forces(Field& section_forces,
                                      Field& contribution_count,
                                      const Field& displacement) {
    (void)section_forces;
    (void)contribution_count;
    (void)displacement;
    return false;
}

} // namespace model
} // namespace fem
