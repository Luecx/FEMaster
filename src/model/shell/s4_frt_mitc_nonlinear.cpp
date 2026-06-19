#include "s4_frt_mitc_nonlinear.h"

#include "../../mattools/vec_util.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace fem::model {

using Vec8 = S4FRTMITC::Vec8;
using Vec24 = S4FRTMITC::Vec24;
using Mat8 = S4FRTMITC::Mat8;
using Mat24 = S4FRTMITC::Mat24;
using Mat8x24 = S4FRTMITC::Mat8x24;
using Mat43 = S4FRTMITC::Mat43;
using Mat42 = S4FRTMITC::Mat42;
using Vec24Mat = S4FRTMITC::Vec24Mat;
using VectorDerivatives = S4FRTMITC::VectorDerivatives;
using ScalarDerivatives = S4FRTMITC::ScalarDerivatives;
using ElementBasis = S4FRTMITC::ElementBasis;
using CurrentState = S4FRTMITC::CurrentState;
using RotationCoefficients = S4FRTMITC::RotationCoefficients;
using StrainData = S4FRTMITC::StrainData;
using NaturalShearData = S4FRTMITC::NaturalShearData;
using EvaluationData = S4FRTMITC::EvaluationData;

using mattools::normalized;
using mattools::skew;

S4FRTMITC::S4FRTMITC(ID id, std::array<ID, 4> nodes)
    : ShellElement<4>(id, nodes)
    , geometry(nodes)
    , integration_scheme_(quadrature::Domain::DOMAIN_ISO_QUAD,
                          quadrature::Order::ORDER_CUBIC) {}


void S4FRTMITC::step_begin() {
    auto data = std::make_unique<ReferenceData>();

    data->X     = init_ref_node_coords();
    data->basis = init_reference_basis(data->X);
    data->area  = Precision(0);

    reference_data_ = std::move(data);
    auto& reference_data = *reference_data_;

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto point = integration_scheme_.get_point(ip);

        reference_data.integration_points[ip] = init_ref_point_data(point.r, point.s, point.w);

        reference_data.area += point.w
                             * reference_data.integration_points[ip].detJ;
    }

    reference_data.tying_points[0] = init_ref_point_data(Precision(0),  Precision(-1), Precision(0));
    reference_data.tying_points[1] = init_ref_point_data(Precision(0),  Precision( 1), Precision(0));
    reference_data.tying_points[2] = init_ref_point_data(Precision(-1), Precision(0),  Precision(0));
    reference_data.tying_points[3] = init_ref_point_data(Precision( 1), Precision(0),  Precision(0));
}

Mat43 S4FRTMITC::init_ref_node_coords() const {
    logging::error(this->_model_data                      != nullptr,
        "S4FRTMITC: no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions_reference != nullptr,
        "S4FRTMITC: POSITION_REFERENCE field is not set");

    const auto& positions = *this->_model_data->positions_reference;
    Mat43 X;
    for (Index i = 0; i < 4; ++i) {
        X.row(i) = positions.row_vec3(this->node_ids[i]).transpose();
    }
    return X;
}

ElementBasis S4FRTMITC::init_reference_basis(const Mat43& X) const {
    ElementBasis basis;

    // Match frt_shell4_analytic.py: construct the element basis from the
    // bilinear mapping derivatives at the element centre.
    const StaticVector<4> dN_dxi  = (StaticVector<4>() <<
        -Precision(0.25),  Precision(0.25), Precision(0.25), -Precision(0.25)
    ).finished();
    const StaticVector<4> dN_deta = (StaticVector<4>() <<
        -Precision(0.25), -Precision(0.25), Precision(0.25),  Precision(0.25)
    ).finished();

    Vec3 X_xi  = Vec3::Zero();
    Vec3 X_eta = Vec3::Zero();
    for (Index i = 0; i < 4; ++i) {
        X_xi  += dN_dxi(i)  * X.row(i).transpose();
        X_eta += dN_deta(i) * X.row(i).transpose();
    }

    basis.e1 = normalized(X_xi);
    basis.e3 = normalized(X_xi.cross(X_eta));
    basis.e2 = normalized(basis.e3.cross(basis.e1));

    for (Index i = 0; i < 4; ++i) {
        basis.d0.row(i) = basis.e3.transpose();
    }

    Vec3 x0 = Vec3::Zero();
    for (Index i = 0; i < 4; ++i) {
        x0 += X.row(i).transpose();
    }
    x0 /= Precision(4);

    for (Index i = 0; i < 4; ++i) {
        const Vec3 dx = X.row(i).transpose() - x0;
        basis.xy(i, 0) = dx.dot(basis.e1);
        basis.xy(i, 1) = dx.dot(basis.e2);
    }

    return basis;
}


void S4FRTMITC::step_end() {
    reference_data_.reset();
}


const S4FRTMITC::ReferenceData& S4FRTMITC::reference_data() const {
    logging::error(
        reference_data_ != nullptr,
        "S4FRTMITC: reference data has not been initialized for element ",
        this->elem_id,
        ". Call step_begin() before using this element path."
    );
    return *reference_data_;
}


const S4FRTMITC::ReferencePointData* S4FRTMITC::cached_reference_point(
    Precision r,
    Precision s) const {
    if (!reference_data_) {
        return nullptr;
    }

    const Precision tol = Precision(1e-14);

    for (const auto& point : reference_data_->integration_points) {
        if (std::abs(point.r - r) <= tol && std::abs(point.s - s) <= tol) {
            return &point;
        }
    }

    for (const auto& point : reference_data_->tying_points) {
        if (std::abs(point.r - r) <= tol && std::abs(point.s - s) <= tol) {
            return &point;
        }
    }

    return nullptr;
}

std::string S4FRTMITC::type_name() const { return "MITC4FRT"; }


std::shared_ptr<SurfaceInterface> S4FRTMITC::surface(int surface_id) {
    return std::make_shared<Surface4>(
        surface_id == 1
            ? std::array<ID, 4>{this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3]}
            : std::array<ID, 4>{this->nodes()[3], this->nodes()[2], this->nodes()[1], this->nodes()[0]}
    );
}


const quadrature::Quadrature& S4FRTMITC::integration_scheme() const {
    return integration_scheme_;
}


StaticVector<4> S4FRTMITC::shape_function(Precision r, Precision s) const {
    StaticVector<4> N;
    N(0) = Precision(0.25) * (Precision(1) - r) * (Precision(1) - s);
    N(1) = Precision(0.25) * (Precision(1) + r) * (Precision(1) - s);
    N(2) = Precision(0.25) * (Precision(1) + r) * (Precision(1) + s);
    N(3) = Precision(0.25) * (Precision(1) - r) * (Precision(1) + s);
    return N;
}


StaticMatrix<4, 2> S4FRTMITC::shape_derivative(Precision r, Precision s) const {
    StaticMatrix<4, 2> dN;
    dN(0, 0) = -Precision(0.25) * (Precision(1) - s);
    dN(1, 0) =  Precision(0.25) * (Precision(1) - s);
    dN(2, 0) =  Precision(0.25) * (Precision(1) + s);
    dN(3, 0) = -Precision(0.25) * (Precision(1) + s);

    dN(0, 1) = -Precision(0.25) * (Precision(1) - r);
    dN(1, 1) = -Precision(0.25) * (Precision(1) + r);
    dN(2, 1) =  Precision(0.25) * (Precision(1) + r);
    dN(3, 1) =  Precision(0.25) * (Precision(1) - r);
    return dN;
}


ShellSection* S4FRTMITC::shell_section() const {
    if (!this->_section) {
        logging::error(false, "Section not set for element ", this->elem_id);
    }
    if (!this->_section->template as<ShellSection>()) {
        logging::error(false, "Section is not a shell section for element ", this->elem_id);
    }
    return this->_section->template as<ShellSection>();
}


// ---------------------------------------------------------------------
// Basic vector/tensor helpers
// ---------------------------------------------------------------------

RotationCoefficients S4FRTMITC::rotation_coefficients(Precision angle_squared) {
    RotationCoefficients coefficients;

    if (angle_squared < Precision(1e-4)) {
        const Precision s2 = angle_squared * angle_squared;
        const Precision s3 = s2 * angle_squared;

        coefficients.a    = Precision(1)
                          - angle_squared / Precision(6)
                          + s2            / Precision(120)
                          - s3            / Precision(5040);
        coefficients.b    = Precision(0.5)
                          - angle_squared / Precision(24)
                          + s2            / Precision(720)
                          - s3            / Precision(40320);
        coefficients.a_s  = -Precision(1) / Precision(6)
                          + angle_squared / Precision(60)
                          - s2            / Precision(1680)
                          + s3            / Precision(90720);
        coefficients.b_s  = -Precision(1) / Precision(24)
                          + angle_squared / Precision(360)
                          - s2            / Precision(13440)
                          + s3            / Precision(907200);
        coefficients.a_ss = Precision(1) / Precision(60)
                          - angle_squared / Precision(840)
                          + s2            / Precision(30240);
        coefficients.b_ss = Precision(1) / Precision(360)
                          - angle_squared / Precision(6720)
                          + s2            / Precision(302400);

        return coefficients;
    }

    const Precision angle       = std::sqrt(angle_squared);
    const Precision angle_cubed = angle_squared * angle;
    const Precision angle_fourth = angle_squared * angle_squared;
    const Precision angle_fifth  = angle_fourth * angle;
    const Precision angle_sixth  = angle_fourth * angle_squared;
    const Precision sin_angle    = std::sin(angle);
    const Precision cos_angle    = std::cos(angle);

    coefficients.a    = sin_angle / angle;
    coefficients.b    = (Precision(1) - cos_angle) / angle_squared;
    coefficients.a_s  = (angle * cos_angle - sin_angle)
                      / (Precision(2) * angle_cubed);
    coefficients.b_s  = (Precision(0.5) * angle * sin_angle
                       + cos_angle - Precision(1))
                      / angle_fourth;
    coefficients.a_ss = (-angle_squared * sin_angle
                       - Precision(3) * angle * cos_angle
                       + Precision(3) * sin_angle)
                      / (Precision(4) * angle_fifth);
    coefficients.b_ss = (angle_squared * cos_angle
                       - Precision(5) * angle * sin_angle
                       - Precision(8) * cos_angle
                       + Precision(8))
                      / (Precision(4) * angle_sixth);

    return coefficients;
}


void S4FRTMITC::rotation_exp_derivatives(
    const Vec3&                                rotation_vector,
    Mat3&                                      rotation,
    std::array<Mat3, 3>&                       first,
    std::array<std::array<Mat3, 3>, 3>&        second) {
    const Precision angle_squared = rotation_vector.squaredNorm();
    const auto      coefficients  = rotation_coefficients(angle_squared);
    const Mat3      omega         = skew(rotation_vector);
    const Mat3      omega_squared = omega * omega;

    const std::array<Mat3, 3> generators{
        skew(Vec3::UnitX()),
        skew(Vec3::UnitY()),
        skew(Vec3::UnitZ())
    };

    std::array<Mat3, 3> omega_derivatives;

    rotation = Mat3::Identity()
             + coefficients.a * omega
             + coefficients.b * omega_squared;

    for (Index i = 0; i < 3; ++i) {
        const Precision a_i = Precision(2)
                            * coefficients.a_s
                            * rotation_vector(i);
        const Precision b_i = Precision(2)
                            * coefficients.b_s
                            * rotation_vector(i);

        omega_derivatives[i] = generators[i] * omega
                             + omega * generators[i];

        first[i] = a_i * omega
                 + coefficients.a * generators[i]
                 + b_i * omega_squared
                 + coefficients.b * omega_derivatives[i];
    }

    for (Index i = 0; i < 3; ++i) {
        const Precision a_i = Precision(2)
                            * coefficients.a_s
                            * rotation_vector(i);
        const Precision b_i = Precision(2)
                            * coefficients.b_s
                            * rotation_vector(i);

        for (Index j = 0; j < 3; ++j) {
            const Precision delta = i == j
                                  ? Precision(1)
                                  : Precision(0);
            const Precision a_j = Precision(2)
                                * coefficients.a_s
                                * rotation_vector(j);
            const Precision b_j = Precision(2)
                                * coefficients.b_s
                                * rotation_vector(j);
            const Precision a_ij = Precision(4)
                                 * coefficients.a_ss
                                 * rotation_vector(i)
                                 * rotation_vector(j)
                                 + Precision(2)
                                 * coefficients.a_s
                                 * delta;
            const Precision b_ij = Precision(4)
                                 * coefficients.b_ss
                                 * rotation_vector(i)
                                 * rotation_vector(j)
                                 + Precision(2)
                                 * coefficients.b_s
                                 * delta;

            second[i][j] = a_ij * omega
                         + a_i * generators[j]
                         + a_j * generators[i]
                         + b_ij * omega_squared
                         + b_i * omega_derivatives[j]
                         + b_j * omega_derivatives[i]
                         + coefficients.b
                         * (generators[i] * generators[j]
                          + generators[j] * generators[i]);
        }
    }
}


Mat3 S4FRTMITC::rotation_exp(const Vec3& rotation_vector) {
    const auto coefficients = rotation_coefficients(
        rotation_vector.squaredNorm()
    );
    const Mat3 omega         = skew(rotation_vector);
    const Mat3 omega_squared = omega * omega;

    return Mat3::Identity()
         + coefficients.a * omega
         + coefficients.b * omega_squared;
}


ScalarDerivatives S4FRTMITC::dot_derivatives(const VectorDerivatives& a,
                                         const VectorDerivatives& b) {
    ScalarDerivatives out;
    out.value = a.value.dot(b.value);
    out.d1    = a.d1 * b.value + b.d1 * a.value;
    out.d2.setZero();

    for (Index i = 0; i < num_dofs; ++i) {
        for (Index j = 0; j < num_dofs; ++j) {
            Precision v = Precision(0);
            v += a.d1.row(i).dot(b.d1.row(j));
            v += a.d1.row(j).dot(b.d1.row(i));
            for (Index c = 0; c < 3; ++c) {
                v += a.d2[c](i, j) * b.value(c);
                v += b.d2[c](i, j) * a.value(c);
            }
            out.d2(i, j) = v;
        }
    }
    return out;
}


ScalarDerivatives S4FRTMITC::scaled(const ScalarDerivatives& a, Precision scale) {
    ScalarDerivatives out;
    out.value = scale * a.value;
    out.d1    = scale * a.d1;
    out.d2    = scale * a.d2;
    return out;
}


VectorDerivatives S4FRTMITC::linear_combination(const std::array<VectorDerivatives, 4>& values,
                                            const StaticVector<4>&                  coeffs) {
    VectorDerivatives out;
    for (Index a = 0; a < 4; ++a) {
        out.value += coeffs(a) * values[a].value;
        out.d1    += coeffs(a) * values[a].d1;
        for (Index c = 0; c < 3; ++c) {
            out.d2[c] += coeffs(a) * values[a].d2[c];
        }
    }
    return out;
}


// ---------------------------------------------------------------------
// Reference/current kinematics
// ---------------------------------------------------------------------

StaticMatrix<4, 6> S4FRTMITC::node_coords_current_6() const {
    logging::error(this->_model_data != nullptr,
                   "S4FRTMITC: no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions != nullptr,
                   "S4FRTMITC: POSITION field is not set");

    const auto& positions = *this->_model_data->positions;
    StaticMatrix<4, 6> q;
    for (Index i = 0; i < 4; ++i) {
        q.row(i) = positions.row_vec6(static_cast<Index>(this->node_ids[i])).transpose();
    }
    return q;
}

S4FRTMITC::ReferencePointData S4FRTMITC::init_ref_point_data(
    Precision           r,
    Precision           s,
    Precision           w) const {
    ReferencePointData point;

    const Mat43&        X     = this->reference_data().X;
    const ElementBasis& basis = this->reference_data().basis;

    point.r     = r;
    point.s     = s;
    point.w     = w;
    point.N     = shape_function(r, s);
    point.dN_rs = shape_derivative(r, s);

    const StaticVector<4> a_nodes = basis.xy.col(0);
    const StaticVector<4> b_nodes = basis.xy.col(1);

    const Precision a_xi  = point.dN_rs.col(0).dot(a_nodes);
    const Precision a_eta = point.dN_rs.col(1).dot(a_nodes);
    const Precision b_xi  = point.dN_rs.col(0).dot(b_nodes);
    const Precision b_eta = point.dN_rs.col(1).dot(b_nodes);

    point.A << a_xi,  b_xi,
               a_eta, b_eta;

    const Precision detA = point.A.determinant();
    logging::error(
        std::abs(detA) > Precision(1e-14),
        "S4FRTMITC: singular reference Jacobian in element ",
        this->elem_id
    );

    point.detJ = std::abs(detA);
    point.invA = point.A.inverse();

    for (Index i = 0; i < num_nodes; ++i) {
        Vec2 rhs;
        rhs << point.dN_rs(i, 0), point.dN_rs(i, 1);
        const Vec2 grad = point.invA * rhs;

        point.dN_da(i) = grad(0);
        point.dN_db(i) = grad(1);
    }

    for (Index i = 0; i < num_nodes; ++i) {
        const Vec3 X_i  = X.row(i).transpose();
        const Vec3 D_i  = basis.d0.row(i).transpose();

        point.X_a   += point.dN_da(i)    * X_i;
        point.X_b   += point.dN_db(i)    * X_i;
        point.X_xi  += point.dN_rs(i, 0) * X_i;
        point.X_eta += point.dN_rs(i, 1) * X_i;
        point.D     += point.N(i)        * D_i;
        point.D_a   += point.dN_da(i)    * D_i;
        point.D_b   += point.dN_db(i)    * D_i;
    }

    return point;
}


CurrentState S4FRTMITC::current_state() const {
    const ElementBasis& basis = reference_data().basis;
    const StaticMatrix<4, 6> positions = node_coords_current_6();

    CurrentState state;
    for (Index i = 0; i < 4; ++i) {
        const Vec3 x               = positions.template block<1, 3>(i, 0).transpose();
        const Vec3 rotation_vector = positions.template block<1, 3>(i, 3).transpose();
        const Vec3 d0              = basis.d0.row(i).transpose();
        const Mat3 rotation        = rotation_exp(rotation_vector);

        state.x.row(i)     = x.transpose();
        state.d.row(i)     = (rotation * d0).transpose();
        state.theta.row(i) = rotation_vector.transpose();
    }
    return state;
}


CurrentState S4FRTMITC::reference_state() const {
    const ElementBasis& basis = reference_data().basis;
    const Mat43 X = reference_data().X;

    CurrentState state;
    state.x     = X;
    state.d     = basis.d0;
    state.theta.setZero();
    return state;
}


CurrentState S4FRTMITC::current_state_from_displacement(
    const Field&        displacement) const {
    logging::error(displacement.domain == FieldDomain::NODE,
                   "S4FRTMITC: displacement field must use NODE domain");
    logging::error(displacement.components >= dofs_per_node,
                   "S4FRTMITC: displacement field requires six components");

    const ElementBasis& basis = reference_data().basis;
    const Mat43 X = reference_data().X;

    CurrentState state;
    for (Index i = 0; i < num_nodes; ++i) {
        const Index node_id = static_cast<Index>(this->node_ids[i]);
        const Vec6  q       = displacement.row_vec6(node_id);
        const Vec3  x       = X.row(i).transpose() + q.template head<3>();
        const Vec3  theta   = q.template tail<3>();
        const Vec3  d0      = basis.d0.row(i).transpose();
        const Mat3  R       = rotation_exp(theta);

        state.x.row(i)     = x.transpose();
        state.d.row(i)     = (R * d0).transpose();
        state.theta.row(i) = theta.transpose();
    }
    return state;
}


Vec24 S4FRTMITC::element_displacement_vector(const Field& displacement) const {
    logging::error(displacement.domain == FieldDomain::NODE,
                   "S4FRTMITC: displacement field must use NODE domain");
    logging::error(displacement.components >= dofs_per_node,
                   "S4FRTMITC: displacement field requires six components");

    Vec24 element_displacement;
    for (Index i = 0; i < num_nodes; ++i) {
        const Index node_id = static_cast<Index>(this->node_ids[i]);
        element_displacement.template segment<6>(dofs_per_node * i) =
            displacement.row_vec6(node_id);
    }
    return element_displacement;
}


Mat3 S4FRTMITC::reference_basis_global() const {
    const ElementBasis& basis = reference_data().basis;
    Mat3 result;
    result.col(0) = basis.e1;
    result.col(1) = basis.e2;
    result.col(2) = basis.e3;
    return result;
}


void S4FRTMITC::shape_gradients_physical(const StaticMatrix<4, 2>& dN_rs,
                              StaticVector<4>&          dN_da,
                              StaticVector<4>&          dN_db,
                              Precision&                detJ,
                              Mat2&                     A) const {
    const ElementBasis& basis = reference_data().basis;
    const StaticVector<4> a_nodes = basis.xy.col(0);
    const StaticVector<4> b_nodes = basis.xy.col(1);

    const Precision a_xi  = dN_rs.col(0).dot(a_nodes);
    const Precision a_eta = dN_rs.col(1).dot(a_nodes);
    const Precision b_xi  = dN_rs.col(0).dot(b_nodes);
    const Precision b_eta = dN_rs.col(1).dot(b_nodes);

    A << a_xi,  b_xi,
         a_eta, b_eta;

    const Precision detA = A.determinant();
    logging::error(std::abs(detA) > Precision(1e-14),
                   "S4FRTMITC: singular reference Jacobian in element ", this->elem_id);

    for (Index i = 0; i < 4; ++i) {
        Vec2 rhs;
        rhs << dN_rs(i, 0), dN_rs(i, 1);
        const Vec2 grad = A.inverse() * rhs;
        dN_da(i) = grad(0);
        dN_db(i) = grad(1);
    }

    detJ = std::abs(detA);
}


std::array<VectorDerivatives, 4> S4FRTMITC::x_derivatives(const CurrentState& state) const {
    std::array<VectorDerivatives, 4> x_nodes;

    for (Index i = 0; i < 4; ++i) {
        const Index k = dofs_per_node * i;
        x_nodes[i].value = state.x.row(i).transpose();
        x_nodes[i].d1(k + 0, 0) = Precision(1);
        x_nodes[i].d1(k + 1, 1) = Precision(1);
        x_nodes[i].d1(k + 2, 2) = Precision(1);
    }
    return x_nodes;
}


std::array<VectorDerivatives, 4> S4FRTMITC::director_derivatives(
    const CurrentState& state) const {
    const ElementBasis& basis = reference_data().basis;
    std::array<VectorDerivatives, 4> directors;

    for (Index i = 0; i < num_nodes; ++i) {
        const Index base            = dofs_per_node * i;
        const Vec3  rotation_vector = state.theta.row(i).transpose();
        const Vec3  d0              = basis.d0.row(i).transpose();

        Mat3 rotation;
        std::array<Mat3, 3> first;
        std::array<std::array<Mat3, 3>, 3> second;

        rotation_exp_derivatives(
            rotation_vector,
            rotation,
            first,
            second
        );

        directors[i].value = rotation * d0;

        for (Index a = 0; a < 3; ++a) {
            const Index ia = base + 3 + a;
            directors[i].d1.row(ia) = (first[a] * d0).transpose();

            for (Index b = 0; b < 3; ++b) {
                const Index ib = base + 3 + b;
                const Vec3  value = second[a][b] * d0;

                for (Index c = 0; c < 3; ++c) {
                    directors[i].d2[c](ia, ib) = value(c);
                }
            }
        }
    }

    return directors;
}


void S4FRTMITC::reference_fields(
    const StaticVector<4>& N,
    const StaticVector<4>& dN_da,
    const StaticVector<4>& dN_db,
    Vec3& X_a,
    Vec3& X_b,
    Vec3& D,
    Vec3& D_a,
    Vec3& D_b
) const {
    const auto& data  = reference_data();
    const auto& X     = data.X;
    const auto& basis = data.basis;

    X_a.setZero();
    X_b.setZero();
    D.setZero();
    D_a.setZero();
    D_b.setZero();

    for (Index i = 0; i < num_nodes; ++i) {
        const Vec3 X_i = X.row(i).transpose();
        const Vec3 D_i = basis.d0.row(i).transpose();

        X_a += dN_da(i) * X_i;
        X_b += dN_db(i) * X_i;
        D   += N(i)      * D_i;
        D_a += dN_da(i) * D_i;
        D_b += dN_db(i) * D_i;
    }
}


EvaluationData S4FRTMITC::create_evaluation_data(
    const CurrentState& state) const {
    EvaluationData data;

    data.state   = state;
    data.x_nodes = x_derivatives(state);
    data.d_nodes = director_derivatives(state);

    const auto& reference = reference_data();
    data.tying.bottom = raw_natural_shear_B_G_at(data, reference.tying_points[0]);
    data.tying.top    = raw_natural_shear_B_G_at(data, reference.tying_points[1]);
    data.tying.left   = raw_natural_shear_B_G_at(data, reference.tying_points[2]);
    data.tying.right  = raw_natural_shear_B_G_at(data, reference.tying_points[3]);

    return data;
}


StrainData S4FRTMITC::raw_strain_B_G_at(
    const EvaluationData&     data,
    const ReferencePointData& point) const {
    StrainData out;

    out.detJ = point.detJ;

    const VectorDerivatives x_a = linear_combination(data.x_nodes, point.dN_da);
    const VectorDerivatives x_b = linear_combination(data.x_nodes, point.dN_db);
    const VectorDerivatives d   = linear_combination(data.d_nodes, point.N);
    const VectorDerivatives d_a = linear_combination(data.d_nodes, point.dN_da);
    const VectorDerivatives d_b = linear_combination(data.d_nodes, point.dN_db);

    std::array<ScalarDerivatives, 8> items;
    items[0] = scaled(dot_derivatives(x_a, x_a), Precision(0.5));
    items[1] = scaled(dot_derivatives(x_b, x_b), Precision(0.5));
    items[2] = dot_derivatives(x_a, x_b);
    items[3] = dot_derivatives(x_a, d_a);
    items[4] = dot_derivatives(x_b, d_b);

    {
        const ScalarDerivatives xa_db = dot_derivatives(x_a, d_b);
        const ScalarDerivatives xb_da = dot_derivatives(x_b, d_a);
        items[5].value = xa_db.value + xb_da.value;
        items[5].d1    = xa_db.d1    + xb_da.d1;
        items[5].d2    = xa_db.d2    + xb_da.d2;
    }

    items[6] = dot_derivatives(x_a, d);
    items[7] = dot_derivatives(x_b, d);

    const std::array<Precision, 8> constants {{
        Precision(0.5) * point.X_a.dot(point.X_a),
        Precision(0.5) * point.X_b.dot(point.X_b),
        point.X_a.dot(point.X_b),
        point.X_a.dot(point.D_a),
        point.X_b.dot(point.D_b),
        point.X_a.dot(point.D_b) + point.X_b.dot(point.D_a),
        point.X_a.dot(point.D),
        point.X_b.dot(point.D)
    }};

    for (Index a = 0; a < num_strains; ++a) {
        out.strain(a) = items[a].value - constants[a];
        out.B.row(a)  = items[a].d1.transpose();
        out.G[a]      = items[a].d2;
    }

    return out;
}


StrainData S4FRTMITC::raw_strain_B_G_at(
    const EvaluationData& data,
    Precision             r,
    Precision             s) const {
    if (const ReferencePointData* point = cached_reference_point(r, s)) {
        return raw_strain_B_G_at(data, *point);
    }

    const ReferencePointData point = init_ref_point_data(r, s, Precision(0));
    return raw_strain_B_G_at(data, point);
}


NaturalShearData S4FRTMITC::raw_natural_shear_B_G_at(
    const EvaluationData&     data,
    const ReferencePointData& point) const {
    NaturalShearData out;

    StaticVector<4> dN_dxi  = point.dN_rs.col(0);
    StaticVector<4> dN_deta = point.dN_rs.col(1);

    const VectorDerivatives x_xi  = linear_combination(data.x_nodes, dN_dxi);
    const VectorDerivatives x_eta = linear_combination(data.x_nodes, dN_deta);
    const VectorDerivatives d     = linear_combination(data.d_nodes, point.N);

    const ScalarDerivatives gamma_xi  = dot_derivatives(x_xi, d);
    const ScalarDerivatives gamma_eta = dot_derivatives(x_eta, d);

    out.shear_nat(0) = gamma_xi.value  - point.X_xi.dot(point.D);
    out.shear_nat(1) = gamma_eta.value - point.X_eta.dot(point.D);

    out.B_nat.row(0) = gamma_xi.d1.transpose();
    out.B_nat.row(1) = gamma_eta.d1.transpose();

    out.G_nat[0] = gamma_xi.d2;
    out.G_nat[1] = gamma_eta.d2;

    return out;
}


NaturalShearData S4FRTMITC::raw_natural_shear_B_G_at(
    const EvaluationData& data,
    Precision             r,
    Precision             s) const {
    if (const ReferencePointData* point = cached_reference_point(r, s)) {
        return raw_natural_shear_B_G_at(data, *point);
    }

    const ReferencePointData point = init_ref_point_data(r, s, Precision(0));
    return raw_natural_shear_B_G_at(data, point);
}


void S4FRTMITC::mitc4_shear_B_G_at(
    const EvaluationData&     data,
    const ReferencePointData& point,
    Vec2&                     shear,
    StaticMatrix<2, num_dofs>& B_shear,
    std::array<Mat24, 2>&     G_shear) const {
    Vec2 g_nat = Vec2::Zero();
    StaticMatrix<2, num_dofs> B_nat = StaticMatrix<2, num_dofs>::Zero();
    std::array<Mat24, 2> G_nat;
    for (auto& matrix : G_nat) {
        matrix.setZero();
    }

    const Precision r = point.r;
    const Precision s = point.s;

    const Precision c_bottom = Precision(0.5) * (Precision(1) - s);
    const Precision c_top    = Precision(0.5) * (Precision(1) + s);
    const Precision c_left   = Precision(0.5) * (Precision(1) - r);
    const Precision c_right  = Precision(0.5) * (Precision(1) + r);

    const auto& bottom = data.tying.bottom;
    const auto& top    = data.tying.top;
    const auto& left   = data.tying.left;
    const auto& right  = data.tying.right;

    g_nat(0)     = c_bottom * bottom.shear_nat(0) + c_top   * top.shear_nat(0);
    B_nat.row(0) = c_bottom * bottom.B_nat.row(0) + c_top   * top.B_nat.row(0);
    G_nat[0]     = c_bottom * bottom.G_nat[0]     + c_top   * top.G_nat[0];

    g_nat(1)     = c_left   * left.shear_nat(1)   + c_right * right.shear_nat(1);
    B_nat.row(1) = c_left   * left.B_nat.row(1)   + c_right * right.B_nat.row(1);
    G_nat[1]     = c_left   * left.G_nat[1]       + c_right * right.G_nat[1];

    shear   = point.invA * g_nat;
    B_shear = point.invA * B_nat;

    for (Index a = 0; a < 2; ++a) {
        G_shear[a].setZero();
        for (Index b = 0; b < 2; ++b) {
            G_shear[a] += point.invA(a, b) * G_nat[b];
        }
    }
}


StrainData S4FRTMITC::strain_B_G_at(
    const EvaluationData&     data,
    const ReferencePointData& point) const {
    StrainData out = raw_strain_B_G_at(data, point);

    Vec2 shear;
    StaticMatrix<2, num_dofs> B_shear;
    std::array<Mat24, 2> G_shear;

    mitc4_shear_B_G_at(data, point, shear, B_shear, G_shear);

    out.strain(6) = shear(0);
    out.strain(7) = shear(1);
    out.B.row(6)  = B_shear.row(0);
    out.B.row(7)  = B_shear.row(1);
    out.G[6]      = G_shear[0];
    out.G[7]      = G_shear[1];

    return out;
}


StrainData S4FRTMITC::strain_B_G_at(
    const EvaluationData& data,
    Precision             r,
    Precision             s) const {
    if (const ReferencePointData* point = cached_reference_point(r, s)) {
        return strain_B_G_at(data, *point);
    }

    const ReferencePointData point = init_ref_point_data(r, s, Precision(0));
    return strain_B_G_at(data, point);
}


StrainData S4FRTMITC::strain_B_G_at(
    const CurrentState& state,
    Precision           r,
    Precision           s) const {
    const EvaluationData data = create_evaluation_data(state);
    return strain_B_G_at(data, r, s);
}


// ---------------------------------------------------------------------
// Material/resultant matrices and EAS
// ---------------------------------------------------------------------

Mat8 S4FRTMITC::resultant_stiffness() const {
    ShellSection* section = shell_section();

    Mat8 H = Mat8::Zero();
    H.template block<6, 6>(0, 0) = section->get_abd();
    H.template block<2, 2>(6, 6) = section->get_shear();

    Precision topo_scale = Precision(1);
    if (this->_model_data && this->_model_data->element_stiffness_scale) {
        auto scale_field = this->_model_data->element_stiffness_scale;
        logging::error(scale_field->components == 1,
                       "Field '", scale_field->name, "': element stiffness scale requires 1 component");
        topo_scale = (*scale_field)(static_cast<Index>(this->elem_id));
    }
    return topo_scale * H;
}


StaticMatrix<S4FRTMITC::num_strains, S4FRTMITC::eas_parameters> S4FRTMITC::eas_matrix(
    Precision r,
    Precision s) const {
    StaticMatrix<num_strains, eas_parameters> M;
    M.setZero();

    if constexpr (eas_parameters >= 4) {
        M(0, 0) = r;
        M(1, 1) = s;
        M(2, 2) = r;
        M(2, 3) = s;
    }
    if constexpr (eas_parameters >= 5) {
        M(0, 4) = r * s;
        M(1, 4) = -r * s;
    }
    return M;
}


StaticVector<S4FRTMITC::eas_parameters> S4FRTMITC::compute_eas_alpha(const CurrentState& state,
                                               const Mat8&         H) const {
    StaticVector<eas_parameters> alpha;
    alpha.setZero();

    if constexpr (eas_parameters == 0) {
        return alpha;
    } else {
        StaticMatrix<eas_parameters, eas_parameters> Kaa;
        StaticVector<eas_parameters> ba;
        Kaa.setZero();
        ba.setZero();

        const EvaluationData data = create_evaluation_data(state);

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto&     point = reference_data().integration_points[ip];
            const StrainData sd   = strain_B_G_at(data, point);
            const Precision  fac  = point.w * sd.detJ;
            const auto       M    = eas_matrix(point.r, point.s);
            const StaticMatrix<eas_parameters, num_strains> MTH = M.transpose() * H;

            Kaa += fac * (MTH * M);
            ba  += fac * (MTH * sd.strain);
        }

        alpha = -Kaa.ldlt().solve(ba);
        return alpha;
    }
}


void S4FRTMITC::material_and_geometric_stiffness(const CurrentState& state,
                                        const Mat8&         H,
                                        Mat24&              Kmat,
                                        Mat24&              Kgeo,
                                        Vec24*              force) const {
    Kmat.setZero();
    Kgeo.setZero();
    if (force) {
        force->setZero();
    }

    if constexpr (eas_parameters > 0) {
        StaticMatrix<eas_parameters, eas_parameters> Kaa;
        StaticVector<eas_parameters>                 ba;
        StaticMatrix<eas_parameters, num_dofs>     Jb;
        std::array<Mat24, eas_parameters>            Hb;
        Kaa.setZero();
        ba.setZero();
        Jb.setZero();
        for (auto& h : Hb) {
            h.setZero();
        }

        const EvaluationData data = create_evaluation_data(state);

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto&      point = reference_data().integration_points[ip];
            const StrainData sd    = strain_B_G_at(data, point);
            const Precision  fac   = point.w * sd.detJ;
            const Vec8       res   = H * sd.strain;

            Kmat += fac * (sd.B.transpose() * H * sd.B);
            for (Index a = 0; a < num_strains; ++a) {
                Kgeo += fac * res(a) * sd.G[a];
            }
            if (force) {
                *force += fac * (sd.B.transpose() * res);
            }

            const auto M = eas_matrix(point.r, point.s);
            const StaticMatrix<eas_parameters, num_strains> MTH = M.transpose() * H;
            Kaa += fac * (MTH * M);
            ba  += fac * (MTH * sd.strain);
            Jb  += fac * (MTH * sd.B);

            for (Index k = 0; k < eas_parameters; ++k) {
                for (Index a = 0; a < num_strains; ++a) {
                    Hb[k] += fac * MTH(k, a) * sd.G[a];
                }
            }
        }

        const auto solver = Kaa.ldlt();
        const StaticVector<eas_parameters> alpha = -solver.solve(ba);

        Kmat -= Jb.transpose() * solver.solve(Jb);
        for (Index k = 0; k < eas_parameters; ++k) {
            Kgeo += alpha(k) * Hb[k];
        }
        if (force) {
            *force += Jb.transpose() * alpha;
        }
    } else {
        const EvaluationData data = create_evaluation_data(state);

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto&     point = reference_data().integration_points[ip];
            const StrainData sd   = strain_B_G_at(data, point);
            const Precision  fac  = point.w * sd.detJ;
            const Vec8       s    = H * sd.strain;

            Kmat += fac * (sd.B.transpose() * H * sd.B);
            for (Index a = 0; a < num_strains; ++a) {
                Kgeo += fac * s(a) * sd.G[a];
            }
            if (force) {
                *force += fac * (sd.B.transpose() * s);
            }
        }
    }
}


Vec8 S4FRTMITC::shell_resultant_at(const CurrentState& state,
                        const Mat8&         H,
                        Precision           r,
                        Precision           s,
                        Vec8*               strain_out) const {
    const StrainData sd = strain_B_G_at(state, r, s);
    Vec8 eps = sd.strain;

    if constexpr (eas_parameters > 0) {
        const auto alpha = compute_eas_alpha(state, H);
        eps += eas_matrix(r, s) * alpha;
    }

    if (strain_out) {
        *strain_out = eps;
    }
    return H * eps;
}


StaticVector<S4FRTMITC::eas_parameters> S4FRTMITC::compute_eas_alpha_linear(
    const Field&        displacement,
    const Mat8&         H) const {
    StaticVector<eas_parameters> alpha;
    alpha.setZero();

    if constexpr (eas_parameters == 0) {
        return alpha;
    } else {
        const CurrentState state0 = reference_state();
        const Vec24        q       = element_displacement_vector(displacement);

        StaticMatrix<eas_parameters, eas_parameters> Kaa;
        StaticVector<eas_parameters>                 ba;
        Kaa.setZero();
        ba.setZero();

        const EvaluationData data0 = create_evaluation_data(state0);

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto&      point = reference_data().integration_points[ip];
            const StrainData sd    = strain_B_G_at(data0, point);
            const Precision  dA    = point.w * sd.detJ;
            const Vec8       eps   = sd.B * q;
            const auto       M     = eas_matrix(point.r, point.s);
            const StaticMatrix<eas_parameters, num_strains> MTH = M.transpose() * H;

            Kaa += dA * (MTH * M);
            ba  += dA * (MTH * eps);
        }

        alpha = -Kaa.ldlt().solve(ba);
        return alpha;
    }
}


Vec8 S4FRTMITC::generalized_strain_at(const Field&        displacement,
                           Precision           r,
                           Precision           s,
                           bool                nonlinear) const {
    const Mat8 H = resultant_stiffness();

    if (nonlinear) {
        const CurrentState state = current_state_from_displacement(displacement);
        const StrainData   data  = strain_B_G_at(state, r, s);
        Vec8               eps   = data.strain;

        if constexpr (eas_parameters > 0) {
            eps += eas_matrix(r, s) * compute_eas_alpha(state, H);
        }
        return eps;
    }

    const CurrentState state0 = reference_state();
    const Vec24        q      = element_displacement_vector(displacement);
    const StrainData   data0  = strain_B_G_at(state0, r, s);
    Vec8               eps    = data0.B * q;

    if constexpr (eas_parameters > 0) {
        eps += eas_matrix(r, s) * compute_eas_alpha_linear(displacement, H);
    }
    return eps;
}


Vec8 S4FRTMITC::generalized_resultant_at(const Field&        displacement,
                              Precision           r,
                              Precision           s,
                              bool                nonlinear,
                              Vec8*               strain_out) const {
    const Vec8 eps = generalized_strain_at(displacement, r, s, nonlinear);

    if (strain_out) {
        *strain_out = eps;
    }
    return resultant_stiffness() * eps;
}


Mat3 S4FRTMITC::current_output_basis(const CurrentState& state,
                          Precision           r,
                          Precision           s) const {
    const StaticVector<4>   N     = shape_function(r, s);
    const StaticMatrix<4,2> dN_rs = shape_derivative(r, s);

    StaticVector<4> dN_da;
    StaticVector<4> dN_db;
    Precision       detJ;
    Mat2            A;
    shape_gradients_physical(dN_rs, dN_da, dN_db, detJ, A);

    Vec3 x_a      = Vec3::Zero();
    Vec3 x_b      = Vec3::Zero();
    Vec3 director = Vec3::Zero();

    for (Index i = 0; i < num_nodes; ++i) {
        x_a      += dN_da(i) * state.x.row(i).transpose();
        x_b      += dN_db(i) * state.x.row(i).transpose();
        director += N(i)     * state.d.row(i).transpose();
    }

    Vec3 e3 = normalized(director);
    Vec3 e1 = x_a - e3 * x_a.dot(e3);

    if (e1.squaredNorm() <= Precision(1e-24)) {
        e1 = x_b - e3 * x_b.dot(e3);
    }

    e1 = normalized(e1);
    Vec3 e2 = normalized(e3.cross(e1));
    e1 = normalized(e2.cross(e3));

    Mat3 output_basis;
    output_basis.col(0) = e1;
    output_basis.col(1) = e2;
    output_basis.col(2) = e3;
    return output_basis;
}


Mat3 S4FRTMITC::deformation_gradient_at(const CurrentState& state,
                             Precision           r,
                             Precision           s,
                             Precision           z) const {
    const StaticVector<4>   N     = shape_function(r, s);
    const StaticMatrix<4,2> dN_rs = shape_derivative(r, s);

    StaticVector<4> dN_da;
    StaticVector<4> dN_db;
    Precision       detJ;
    Mat2            A;
    shape_gradients_physical(dN_rs, dN_da, dN_db, detJ, A);

    Vec3 X_a;
    Vec3 X_b;
    Vec3 D;
    Vec3 D_a;
    Vec3 D_b;
    reference_fields(N, dN_da, dN_db, X_a, X_b, D, D_a, D_b);

    Vec3 x_a = Vec3::Zero();
    Vec3 x_b = Vec3::Zero();
    Vec3 d   = Vec3::Zero();
    Vec3 d_a = Vec3::Zero();
    Vec3 d_b = Vec3::Zero();

    for (Index i = 0; i < num_nodes; ++i) {
        const Vec3 x_i = state.x.row(i).transpose();
        const Vec3 d_i = state.d.row(i).transpose();

        x_a += dN_da(i) * x_i;
        x_b += dN_db(i) * x_i;
        d   += N(i)     * d_i;
        d_a += dN_da(i) * d_i;
        d_b += dN_db(i) * d_i;
    }

    Mat3 reference_covariant;
    reference_covariant.col(0) = X_a + z * D_a;
    reference_covariant.col(1) = X_b + z * D_b;
    reference_covariant.col(2) = D;

    Mat3 current_covariant;
    current_covariant.col(0) = x_a + z * d_a;
    current_covariant.col(1) = x_b + z * d_b;
    current_covariant.col(2) = d;

    const Precision reference_det = reference_covariant.determinant();
    logging::error(std::abs(reference_det) > Precision(1e-14),
                   "S4FRTMITC: singular reference shell basis during stress recovery");

    return current_covariant * reference_covariant.inverse();
}


void S4FRTMITC::physical_stress_strain_at(const Field&        displacement,
                               Precision           r,
                               Precision           s,
                               Precision           zeta,
                               bool                nonlinear,
                               Vec6&               strain_out,
                               Vec6&               stress_out) const {
    Vec8 generalized_strain;
    const Vec8 resultants = generalized_resultant_at(
        displacement,
        r,
        s,
        nonlinear,
        &generalized_strain
    );

    const Precision h = this->get_section()->thickness_;
    const Precision z = Precision(0.5) * h * zeta;

    const Vec3 plane_strain =
        generalized_strain.template segment<3>(0)
        + z * generalized_strain.template segment<3>(3);

    const Vec3 plane_stress =
        resultants.template segment<3>(0) / h
        + z * (Precision(12) / (h * h * h))
        * resultants.template segment<3>(3);

    const Vec2 shear_strain = generalized_strain.template segment<2>(6);
    const Vec2 shear_stress = resultants.template segment<2>(6) / h;

    Vec6 strain_local = Vec6::Zero();
    Vec6 stress_local = Vec6::Zero();

    // Solver Voigt ordering:
    // [11, 22, 33, 23, 13, 12].
    strain_local(0) = plane_strain(0);
    strain_local(1) = plane_strain(1);
    strain_local(3) = shear_strain(1);
    strain_local(4) = shear_strain(0);
    strain_local(5) = plane_strain(2);

    stress_local(0) = plane_stress(0);
    stress_local(1) = plane_stress(1);
    stress_local(3) = shear_stress(1);
    stress_local(4) = shear_stress(0);
    stress_local(5) = plane_stress(2);

    const Mat3 reference_basis_matrix = reference_basis_global();

    const Mat3 green_lagrange_global =
        reference_basis_matrix
        * GreenLagrangeStrain(strain_local).tensor()
        * reference_basis_matrix.transpose();

    const Mat3 second_pk_global =
        reference_basis_matrix
        * PK2Stress(stress_local).tensor()
        * reference_basis_matrix.transpose();

    strain_out = GreenLagrangeStrain(green_lagrange_global).voigt();

    if (!nonlinear) {
        stress_out = PK2Stress(second_pk_global).voigt();
        return;
    }

    const CurrentState state =
        current_state_from_displacement(displacement);
    const Mat3 F = deformation_gradient_at(state, r, s, z);
    const Precision J = F.determinant();

    logging::error(J > Precision(0) && std::isfinite(J),
                   "S4FRTMITC: invalid deformation gradient during stress recovery in element ",
                   this->elem_id,
                   ", J = ", J);

    const Mat3 cauchy_global =
        (F * second_pk_global * F.transpose()) / J;

    stress_out = CauchyStress(cauchy_global).voigt();
}


// ---------------------------------------------------------------------
// Drill stabilization and StructuralElement interface
// ---------------------------------------------------------------------

Precision S4FRTMITC::reference_area() const {
    return reference_data().area;
}


Precision S4FRTMITC::drill_stiffness_per_node(const Mat8& H) const {
    const Precision area = reference_area();
    const Precision shear_scale = Precision(0.5)
                                  * (std::abs(H(6, 6)) + std::abs(H(7, 7)));
    return drill_scale * shear_scale * area / Precision(num_nodes);
}


void S4FRTMITC::add_drill_stiffness(Mat24&              stiffness_matrix,
                         const Mat8&         H) const {
    const ElementBasis& basis = reference_data().basis;
    const Precision kd = drill_stiffness_per_node(H);

    for (Index i = 0; i < num_nodes; ++i) {
        const Vec3 d0    = basis.d0.row(i).transpose();
        const Mat3 block = kd * (d0 * d0.transpose());

        stiffness_matrix.template block<3, 3>(
            dofs_per_node * i + 3,
            dofs_per_node * i + 3
        ) += block;
    }
}


void S4FRTMITC::add_drill_force(Vec24&              internal_force,
                     const CurrentState& state,
                     const Mat8&         H) const {
    const ElementBasis& basis = reference_data().basis;
    const Precision kd = drill_stiffness_per_node(H);

    for (Index i = 0; i < num_nodes; ++i) {
        const Vec3 d0              = basis.d0.row(i).transpose();
        const Vec3 rotation_vector = state.theta.row(i).transpose();
        const Precision drill      = rotation_vector.dot(d0);

        internal_force.template segment<3>(
            dofs_per_node * i + 3
        ) += kd * drill * d0;
    }
}


MapMatrix S4FRTMITC::stiffness(Precision* buffer) {
    const CurrentState state = current_state();
    const Mat8         H     = resultant_stiffness();

    Mat24 material_stiffness;
    Mat24 geometric_dummy;

    material_and_geometric_stiffness(
        state,
        H,
        material_stiffness,
        geometric_dummy,
        nullptr
    );

    add_drill_stiffness(material_stiffness, H);
    material_stiffness = Precision(0.5)
                       * (material_stiffness + material_stiffness.transpose());

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = material_stiffness;
    return mapped;
}


MapMatrix S4FRTMITC::stiffness_geom(Precision*   buffer,
                         const Field& ip_stress,
                         int          ip_start_idx) {
    logging::error(
        ip_stress.components >= num_strains,
        "S4FRTMITC: geometric stiffness requires eight shell resultants "
        "[N11,N22,N12,M11,M22,M12,Q13,Q23]"
    );

    const CurrentState   state = current_state();
    const EvaluationData data  = create_evaluation_data(state);

    Mat24 geometric_stiffness;
    geometric_stiffness.setZero();

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto&      point  = reference_data().integration_points[ip];
        const StrainData strain = strain_B_G_at(data, point);
        const Precision  factor = point.w * strain.detJ;
        const Index      row    = static_cast<Index>(ip_start_idx) + ip;

        Vec8 resultants;
        for (Index a = 0; a < num_strains; ++a) {
            resultants(a) = ip_stress(row, a);
        }

        // Rao's finite-rotation tangent contracts the complete generalized
        // resultant vector with the second strain variation. Consequently,
        // membrane, bending and transverse-shear contributions all enter.
        for (Index a = 0; a < num_strains; ++a) {
            geometric_stiffness += factor * resultants(a) * strain.G[a];
        }
    }

    geometric_stiffness = Precision(0.5)
                        * (geometric_stiffness + geometric_stiffness.transpose());

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = geometric_stiffness;
    return mapped;
}


void S4FRTMITC::compute_internal_force_nonlinear(Field&       node_forces,
                                      const Field& ip_stress,
                                      int          ip_offset) {
    logging::error(
        ip_stress.components >= num_strains,
        "S4FRTMITC: internal force requires eight shell resultants "
        "[N11,N22,N12,M11,M22,M12,Q13,Q23]"
    );

    const CurrentState   state = current_state();
    const EvaluationData data  = create_evaluation_data(state);
    const Mat8           H     = resultant_stiffness();

    Vec24 internal_force = Vec24::Zero();

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto&      point  = reference_data().integration_points[ip];
        const StrainData strain = strain_B_G_at(data, point);
        const Precision  factor = point.w * strain.detJ;
        const Index      row    = static_cast<Index>(ip_offset) + ip;

        Vec8 resultants;
        for (Index a = 0; a < num_strains; ++a) {
            resultants(a) = ip_stress(row, a);
        }

        internal_force += factor * (strain.B.transpose() * resultants);
    }

    add_drill_force(internal_force, state, H);

    for (Index a = 0; a < num_nodes; ++a) {
        const Index node_id = static_cast<Index>(this->node_ids[a]);

        for (Index d = 0; d < dofs_per_node; ++d) {
            node_forces(node_id, d) +=
                internal_force(dofs_per_node * a + d);
        }
    }
}


void S4FRTMITC::compute_stress_strain(Field*           strain,
                           Field*           stress_out,
                           const Field&     displacement,
                           const RowMatrix& rst,
                           int              offset,
                           bool             use_green_lagrange_nl) {
    logging::error(strain != nullptr || stress_out != nullptr,
                   "S4FRTMITC: compute_stress_strain requires at least one output field");
    logging::error(rst.cols() >= 3,
                   "S4FRTMITC: stress/strain coordinates require r,s,t columns");

    for (Index n = 0; n < rst.rows(); ++n) {
        Vec6 strain_value;
        Vec6 stress_value;

        physical_stress_strain_at(
            displacement,
            rst(n, 0),
            rst(n, 1),
            rst(n, 2),
            use_green_lagrange_nl,
            strain_value,
            stress_value
        );

        const Index row = static_cast<Index>(offset) + n;

        if (strain) {
            for (Index component = 0; component < strain->components; ++component) {
                (*strain)(row, component) =
                    component < 6 ? strain_value(component) : Precision(0);
            }
        }

        if (stress_out) {
            for (Index component = 0; component < stress_out->components; ++component) {
                (*stress_out)(row, component) =
                    component < 6 ? stress_value(component) : Precision(0);
            }
        }
    }
}


void S4FRTMITC::compute_stress_state(Field&       stress_state,
                          const Field& displacement,
                          int          offset,
                          bool         use_green_lagrange_nl) {
    logging::error(stress_state.components >= num_strains,
                   "S4FRTMITC: stress state requires at least eight components "
                   "[N11,N22,N12,M11,M22,M12,Q13,Q23]");

    const Mat8 H = resultant_stiffness();

    if (use_green_lagrange_nl) {
        const CurrentState   state = current_state_from_displacement(displacement);
        const EvaluationData data  = create_evaluation_data(state);

        StaticVector<eas_parameters> alpha;
        alpha.setZero();
        if constexpr (eas_parameters > 0) {
            alpha = compute_eas_alpha(state, H);
        }

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto&      point = reference_data().integration_points[ip];
            const StrainData sd    = strain_B_G_at(data, point);
            Vec8             eps   = sd.strain;

            if constexpr (eas_parameters > 0) {
                eps += eas_matrix(point.r, point.s) * alpha;
            }

            const Vec8  values = H * eps;
            const Index row    = static_cast<Index>(offset) + ip;

            for (Index component = 0; component < stress_state.components; ++component) {
                stress_state(row, component) =
                    component < num_strains ? values(component) : Precision(0);
            }
        }

        return;
    }

    const CurrentState   state = reference_state();
    const EvaluationData data  = create_evaluation_data(state);
    const Vec24          q     = element_displacement_vector(displacement);

    StaticVector<eas_parameters> alpha;
    alpha.setZero();
    if constexpr (eas_parameters > 0) {
        alpha = compute_eas_alpha_linear(displacement, H);
    }

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto&      point = reference_data().integration_points[ip];
        const StrainData sd    = strain_B_G_at(data, point);
        Vec8             eps   = sd.B * q;

        if constexpr (eas_parameters > 0) {
            eps += eas_matrix(point.r, point.s) * alpha;
        }

        const Vec8  values = H * eps;
        const Index row    = static_cast<Index>(offset) + ip;

        for (Index component = 0; component < stress_state.components; ++component) {
            stress_state(row, component) =
                component < num_strains ? values(component) : Precision(0);
        }
    }
}

RowMatrix S4FRTMITC::stress_strain_ip_rst() {
    const auto& scheme = integration_scheme();
    RowMatrix rst(scheme.count(), 3);
    rst.setZero();
    for (Index i = 0; i < scheme.count(); ++i) {
        rst(i, 0) = scheme.get_point(i).r;
        rst(i, 1) = scheme.get_point(i).s;
        rst(i, 2) = Precision(0);
    }
    return rst;
}


RowMatrix S4FRTMITC::stress_strain_nodal_rst() {
    RowMatrix rst(4, 3);
    rst << Precision(-1), Precision(-1), Precision(0),
            Precision( 1), Precision(-1), Precision(0),
            Precision( 1), Precision( 1), Precision(0),
            Precision(-1), Precision( 1), Precision(0);
    return rst;
}


// ---------------------------------------------------------------------
// Mass, volume and load integration helpers
// ---------------------------------------------------------------------

Precision S4FRTMITC::volume() {
    return this->get_section()->thickness_ * reference_area();
}


StaticMatrix<4, 4> S4FRTMITC::integrate_NNt() const {
    StaticMatrix<4, 4> M;
    M.setZero();

    for (const auto& point : reference_data().integration_points) {
        M += point.w * point.detJ * (point.N * point.N.transpose());
    }
    return M;
}


MapMatrix S4FRTMITC::mass(Precision* buffer) {
    const ElementBasis& basis = reference_data().basis;
    const Precision    rho   = this->get_density(false);

    Mat24 mass_matrix;
    mass_matrix.setZero();

    if (rho > Precision(0)) {
        const Precision h            = this->get_section()->thickness_;
        const Precision mass_surface = rho * h;
        const Precision inertia_rot  =
            rho * h * h * h / Precision(12);
        const auto MNN = integrate_NNt();

        const Vec3 d0 = basis.e3;
        const Mat3 tangent_projector =
            Mat3::Identity() - d0 * d0.transpose();
        const Mat3 drill_projector = d0 * d0.transpose();

        for (Index i = 0; i < num_nodes; ++i) {
            for (Index j = 0; j < num_nodes; ++j) {
                const Precision mij_t = mass_surface * MNN(i, j);
                const Precision mij_r = inertia_rot  * MNN(i, j);

                mass_matrix.template block<3, 3>(
                    dofs_per_node * i,
                    dofs_per_node * j
                ) += mij_t * Mat3::Identity();

                mass_matrix.template block<3, 3>(
                    dofs_per_node * i + 3,
                    dofs_per_node * j + 3
                ) += mij_r * tangent_projector;
            }
        }

        Precision average_nodal_mass = Precision(0);

        for (Index i = 0; i < num_nodes; ++i) {
            average_nodal_mass += mass_surface * MNN(i, i);
        }

        average_nodal_mass /= Precision(num_nodes);

        const Precision drill_inertia =
            Precision(1e-6) * average_nodal_mass;

        for (Index i = 0; i < num_nodes; ++i) {
            mass_matrix.template block<3, 3>(
                dofs_per_node * i + 3,
                dofs_per_node * i + 3
            ) += drill_inertia * drill_projector;
        }
    }

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = mass_matrix;
    return mapped;
}


Vec3 S4FRTMITC::global_point_current(Precision r, Precision s) const {
    const auto N   = shape_function(r, s);
    const auto pos = node_coords_current_6();
    Vec3 x = Vec3::Zero();
    for (Index i = 0; i < 4; ++i) {
        x += N(i) * pos.template block<1, 3>(i, 0).transpose();
    }
    return x;
}


Precision S4FRTMITC::integrate_scalar_field(bool scale_by_density, const ScalarField& field) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

    Precision result = Precision(0);
    for (const auto& point : reference_data().integration_points) {
        result += field(global_point_current(point.r, point.s))
                * (rho * h * point.w * point.detJ);
    }
    return result;
}


Vec3 S4FRTMITC::integrate_vector_field(bool scale_by_density, const VecField& field) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

    Vec3 result = Vec3::Zero();
    for (const auto& point : reference_data().integration_points) {
        result += field(global_point_current(point.r, point.s))
                * (rho * h * point.w * point.detJ);
    }
    return result;
}


void S4FRTMITC::integrate_vector_field(Field& node_loads, bool scale_by_density, const VecField& field) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

    for (const auto& point : reference_data().integration_points) {
        const Vec3 f = field(global_point_current(point.r, point.s))
                     * (rho * h * point.w * point.detJ);

        for (Index i = 0; i < num_nodes; ++i) {
            const Index node_id = static_cast<Index>(this->node_ids[i]);
            for (Index d = 0; d < 3; ++d) {
                node_loads(node_id, d) += point.N(i) * f(d);
            }
        }
    }
}


Mat3 S4FRTMITC::integrate_tensor_field(bool scale_by_density, const TenField& field) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

    Mat3 result = Mat3::Zero();
    for (const auto& point : reference_data().integration_points) {
        result += field(global_point_current(point.r, point.s))
                * (rho * h * point.w * point.detJ);
    }
    return result;
}


void S4FRTMITC::compute_compliance(Field& displacement, Field& result) {
    Precision buffer[num_dofs * num_dofs];
    MapMatrix K = stiffness(buffer);
    const Vec24     u = element_displacement_vector(displacement);

    result(static_cast<Index>(this->elem_id), 0) = u.dot(K * u);
}


void S4FRTMITC::compute_compliance_angle_derivative(Field& displacement,
                                         Field& result) {
    (void) displacement;
    (void) result;
}


bool S4FRTMITC::compute_shear_flow(Field&       shear_flow,
                        const Field& displacement,
                        int          offset) {
    (void) shear_flow;
    (void) displacement;
    (void) offset;
    return false;
}


bool S4FRTMITC::compute_beam_section_forces(Field&       section_forces,
                                 const Field& displacement,
                                 int          offset) {
    (void) section_forces;
    (void) displacement;
    (void) offset;
    return false;
}


bool S4FRTMITC::compute_shell_section_forces(Field&       resultants,
                                  Field&       contribution_count,
                                  const Field& displacement) {
    logging::error(resultants.components >= num_strains,
                   "S4FRTMITC: shell section forces require eight components "
                   "[N11,N22,N12,M11,M22,M12,Q13,Q23]");

    const RowMatrix rst = stress_strain_nodal_rst();

    for (Index i = 0; i < num_nodes; ++i) {
        const Vec8 values = generalized_resultant_at(
            displacement,
            rst(i, 0),
            rst(i, 1),
            true
        );
        const Index node_id = static_cast<Index>(this->node_ids[i]);

        for (Index component = 0; component < num_strains; ++component) {
            resultants(node_id, component) += values(component);
        }
        contribution_count(node_id, 0) += Precision(1);
    }
    return true;
}

} // namespace fem::model
