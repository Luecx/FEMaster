#include "qspt.h"

namespace fem::model {
namespace {
constexpr Precision kGeomTol = Precision(1e-12);

inline Vec3 row_as_vec3(const QSPT::NodeCoords& coords, Index i) {
    return coords.row(i).transpose();
}

inline Precision leverage_to_edge(const Vec3& origin,
                                  const Vec3& midpoint,
                                  const Vec3& edge,
                                  Precision   edge_length) {
    const Vec3 r = midpoint - origin;
    return r.cross(edge).norm() / edge_length;
}
} // namespace

QSPT::QSPT(ID p_elem_id, std::array<ID, 4> p_node_ids)
    : ShellElement<4>(p_elem_id, p_node_ids)
    , geometry(p_node_ids)
    , integration_scheme_(quadrature::DOMAIN_ISO_QUAD, quadrature::ORDER_CUBIC) {}

SurfacePtr QSPT::surface(ID surface_id) {
    return std::make_shared<Surface4>(
        surface_id == 1
            ? std::array<ID, 4>{this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3]}
            : std::array<ID, 4>{this->nodes()[3], this->nodes()[2], this->nodes()[1], this->nodes()[0]});
}

const quadrature::Quadrature& QSPT::integration_scheme() const {
    return integration_scheme_;
}

Precision QSPT::volume() {
    logging::error(this->_model_data != nullptr, "no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions != nullptr, "positions field not set in model data");
    return this->get_section()->thickness_ * geometry.area(*this->_model_data->positions);
}

QSPT::GeometryData QSPT::geometry_data() {
    GeometryData data;
    data.coords = this->node_coords_global();

    data.midpoints[0] = Precision(0.5) * (row_as_vec3(data.coords, 0) + row_as_vec3(data.coords, 1));
    data.midpoints[1] = Precision(0.5) * (row_as_vec3(data.coords, 1) + row_as_vec3(data.coords, 2));
    data.midpoints[2] = Precision(0.5) * (row_as_vec3(data.coords, 2) + row_as_vec3(data.coords, 3));
    data.midpoints[3] = Precision(0.5) * (row_as_vec3(data.coords, 3) + row_as_vec3(data.coords, 0));

    data.edges[0] = row_as_vec3(data.coords, 1) - row_as_vec3(data.coords, 0);
    data.edges[1] = row_as_vec3(data.coords, 2) - row_as_vec3(data.coords, 1);
    data.edges[2] = row_as_vec3(data.coords, 3) - row_as_vec3(data.coords, 2);
    data.edges[3] = row_as_vec3(data.coords, 0) - row_as_vec3(data.coords, 3);

    for (Index i = 0; i < 4; ++i) {
        data.edge_lengths[i] = data.edges[i].norm();
        logging::error(data.edge_lengths[i] > kGeomTol,
                       "QSPT edge length is zero for element ", this->elem_id, " at edge ", i);
    }

    const auto node1_lev = [&]() {
        std::array<Precision, 4> values {};
        const Vec3 origin = row_as_vec3(data.coords, 0);
        for (Index i = 0; i < 4; ++i) {
            values[i] = leverage_to_edge(origin, data.midpoints[i], data.edges[i], data.edge_lengths[i]);
        }
        return values;
    }();
    const auto node2_lev = [&]() {
        std::array<Precision, 4> values {};
        const Vec3 origin = row_as_vec3(data.coords, 1);
        for (Index i = 0; i < 4; ++i) {
            values[i] = leverage_to_edge(origin, data.midpoints[i], data.edges[i], data.edge_lengths[i]);
        }
        return values;
    }();
    const auto node3_lev = [&]() {
        std::array<Precision, 4> values {};
        const Vec3 origin = row_as_vec3(data.coords, 2);
        for (Index i = 0; i < 4; ++i) {
            values[i] = leverage_to_edge(origin, data.midpoints[i], data.edges[i], data.edge_lengths[i]);
        }
        return values;
    }();

    data.h12 = node1_lev[1];
    data.h13 = node1_lev[2];
    data.h23 = node2_lev[2];
    data.h24 = node2_lev[3];
    data.h34 = node3_lev[3];
    data.h31 = node3_lev[0];

    logging::error(data.h12 > kGeomTol && data.h13 > kGeomTol &&
                   data.h23 > kGeomTol && data.h24 > kGeomTol &&
                   data.h34 > kGeomTol && data.h31 > kGeomTol,
                   "QSPT leverage became singular for element ", this->elem_id);

    const Precision a1 = data.edges[0].cross(-data.edges[3]).norm();
    const Precision a2 = (-data.edges[1]).cross(data.edges[2]).norm();
    data.area = Precision(0.5) * (a1 + a2);

    logging::error(data.area > kGeomTol, "QSPT area is zero for element ", this->elem_id);

    return data;
}

Precision QSPT::effective_density() {
    return this->get_density(false);
}

Precision QSPT::effective_shear_modulus() {
    const auto memb = this->get_elasticity()->get_2d();
    const Precision g = memb(2, 2);
    logging::error(std::abs(g) > kGeomTol,
                   "QSPT requires a non-zero in-plane shear modulus for element ", this->elem_id);
    return g;
}

QSPT::ShearState QSPT::shear_state() {
    const GeometryData data = geometry_data();
    const Precision h = data.h12 * data.h23 / (data.h13 * data.h24);
    logging::error(h > kGeomTol, "QSPT metric parameter h became invalid for element ", this->elem_id);

    ShearState state;
    state.edge_lengths = data.edge_lengths;

    const Precision base = std::sqrt(data.edge_lengths[1] * data.edge_lengths[3] / h);
    state.q(0) = -data.h34 / data.h31 * h * base;
    state.q(1) = base;
    state.q(2) = -data.h12 / data.h13 * base;
    state.q(3) = std::sqrt(data.edge_lengths[1] * data.edge_lengths[3] * h);

    for (Index i = 0; i < 4; ++i) {
        const Index prev = (i + 3) % 4;
        const Vec3  term =
            Precision(0.5) * (data.edges[prev] / data.edge_lengths[prev] * state.q(prev) +
                              data.edges[i] / data.edge_lengths[i] * state.q(i));
        state.fn.template segment<3>(3 * i) = term;
    }

    state.flexibility = data.area / (effective_shear_modulus() * this->get_section()->thickness_);
    logging::error(state.flexibility > kGeomTol,
                   "QSPT flexibility became singular for element ", this->elem_id);
    return state;
}

QSPT::StiffnessMatrix QSPT::stiffness_impl() {
    const ShearState state = shear_state();
    return (Precision(1) / state.flexibility) * (state.fn * state.fn.transpose());
}

QSPT::MassMatrix QSPT::mass_impl() {
    MassMatrix M = MassMatrix::Zero();
    const Precision rho = effective_density();
    if (rho <= Precision(0)) {
        return M;
    }

    const Precision thickness = this->get_section()->thickness_;
    const NodeCoords coords = this->node_coords_global();

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto p = integration_scheme_.get_point(ip);
        const ShapeFunction N = geometry.shape_function(p.r, p.s);
        const auto J = geometry.jacobian(coords, p.r, p.s);
        const Precision dA = J.col(0).cross(J.col(1)).norm();
        const Precision scale = rho * thickness * p.w * dA;

        for (Index i = 0; i < 4; ++i) {
            for (Index j = 0; j < 4; ++j) {
                const Precision mij = scale * N(i) * N(j);
                for (Index d = 0; d < 3; ++d) {
                    M(3 * i + d, 3 * j + d) += mij;
                }
            }
        }
    }

    return M;
}

MapMatrix QSPT::stiffness(Precision* buffer) {
    MapMatrix mapped(buffer, 12, 12);
    mapped = stiffness_impl();
    return mapped;
}

MapMatrix QSPT::stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) {
    (void) ip_stress;
    (void) ip_start_idx;
    MapMatrix mapped(buffer, 12, 12);
    mapped.setZero();
    return mapped;
}

MapMatrix QSPT::mass(Precision* buffer) {
    MapMatrix mapped(buffer, 12, 12);
    mapped = mass_impl();
    return mapped;
}

void QSPT::compute_stress_strain(Field& ip_stress,
                                 Field& ip_strain,
                                 Field& displacement,
                                 int ip_offset) {
    (void) ip_stress;
    (void) ip_strain;
    (void) displacement;
    (void) ip_offset;
}

void QSPT::compute_stress_strain_nodal(Field& displacement,
                                       Field& stress,
                                       Field& strain) {
    (void) displacement;
    (void) stress;
    (void) strain;
}

Stresses QSPT::stress(Field& displacement, std::vector<Vec3>& rst) {
    (void) displacement;
    (void) rst;
    return {};
}

Strains QSPT::strain(Field& displacement, std::vector<Vec3>& rst) {
    (void) displacement;
    (void) rst;
    return {};
}

StaticVector<12> QSPT::displacement_vector(Field& displacement) {
    const auto u_mat = this->template nodal_data<3>(displacement);
    StaticVector<12> u = StaticVector<12>::Zero();
    for (Index i = 0; i < 4; ++i) {
        u.template segment<3>(3 * i) = u_mat.row(i).transpose();
    }
    return u;
}

std::vector<Precision> QSPT::shear_flow(Field& displacement) {
    const ShearState state = shear_state();
    const StaticVector<12> u = displacement_vector(displacement);
    const Precision beta = (state.fn.dot(u)) / state.flexibility;

    std::vector<Precision> result(4, Precision(0));
    for (Index i = 0; i < 4; ++i) {
        result[static_cast<std::size_t>(i)] = state.q(i) / state.edge_lengths[i] * beta;
    }
    return result;
}

Precision QSPT::integrate_scalar_field(bool scale_by_density, const ScalarField& field) {
    const NodeCoords coords = this->node_coords_global();
    const Precision density_scale = scale_by_density ? effective_density() : Precision(1);
    const Precision thickness = this->get_section()->thickness_;

    Precision result = Precision(0);
    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto p = integration_scheme_.get_point(ip);
        const ShapeFunction N = geometry.shape_function(p.r, p.s);
        const auto J = geometry.jacobian(coords, p.r, p.s);
        const Precision dA = J.col(0).cross(J.col(1)).norm();

        Vec3 x = Vec3::Zero();
        for (Index i = 0; i < 4; ++i) {
            x += N(i) * row_as_vec3(coords, i);
        }

        result += field(x) * density_scale * thickness * p.w * dA;
    }
    return result;
}

Vec3 QSPT::integrate_vector_field(bool scale_by_density, const VecField& field) {
    const NodeCoords coords = this->node_coords_global();
    const Precision density_scale = scale_by_density ? effective_density() : Precision(1);
    const Precision thickness = this->get_section()->thickness_;

    Vec3 result = Vec3::Zero();
    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto p = integration_scheme_.get_point(ip);
        const ShapeFunction N = geometry.shape_function(p.r, p.s);
        const auto J = geometry.jacobian(coords, p.r, p.s);
        const Precision dA = J.col(0).cross(J.col(1)).norm();

        Vec3 x = Vec3::Zero();
        for (Index i = 0; i < 4; ++i) {
            x += N(i) * row_as_vec3(coords, i);
        }

        result += field(x) * density_scale * thickness * p.w * dA;
    }
    return result;
}

Mat3 QSPT::integrate_tensor_field(bool scale_by_density, const TenField& field) {
    const NodeCoords coords = this->node_coords_global();
    const Precision density_scale = scale_by_density ? effective_density() : Precision(1);
    const Precision thickness = this->get_section()->thickness_;

    Mat3 result = Mat3::Zero();
    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto p = integration_scheme_.get_point(ip);
        const ShapeFunction N = geometry.shape_function(p.r, p.s);
        const auto J = geometry.jacobian(coords, p.r, p.s);
        const Precision dA = J.col(0).cross(J.col(1)).norm();

        Vec3 x = Vec3::Zero();
        for (Index i = 0; i < 4; ++i) {
            x += N(i) * row_as_vec3(coords, i);
        }

        result += field(x) * density_scale * thickness * p.w * dA;
    }
    return result;
}

void QSPT::integrate_vec_field(Field& node_loads,
                               bool scale_by_density,
                               const VecField& field) {
    const NodeCoords coords = this->node_coords_global();
    const Precision density_scale = scale_by_density ? effective_density() : Precision(1);
    const Precision thickness = this->get_section()->thickness_;

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto p = integration_scheme_.get_point(ip);
        const ShapeFunction N = geometry.shape_function(p.r, p.s);
        const auto J = geometry.jacobian(coords, p.r, p.s);
        const Precision dA = J.col(0).cross(J.col(1)).norm();

        Vec3 x = Vec3::Zero();
        for (Index i = 0; i < 4; ++i) {
            x += N(i) * row_as_vec3(coords, i);
        }

        const Vec3 value = field(x) * density_scale * thickness * p.w * dA;
        for (Index i = 0; i < 4; ++i) {
            const ID node_id = this->node_ids[i];
            node_loads(node_id, 0) += N(i) * value(0);
            node_loads(node_id, 1) += N(i) * value(1);
            node_loads(node_id, 2) += N(i) * value(2);
        }
    }
}

void QSPT::compute_compliance(Field& displacement, Field& result) {
    const StaticVector<12> u = displacement_vector(displacement);
    const Precision c = u.dot(stiffness_impl() * u);
    result(this->elem_id, 0) = c;
}

void QSPT::compute_compliance_angle_derivative(Field& displacement, Field& result) {
    (void) displacement;
    (void) result;
}
} // namespace fem::model
