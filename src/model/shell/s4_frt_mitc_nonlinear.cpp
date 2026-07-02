#include "s4_frt_mitc_nonlinear.h"

#include "../../math/so3.h"
#include "../../math/vec_util.h"

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
using CurrentState = S4FRTMITC::CurrentState;
using EvaluationData = S4FRTMITC::EvaluationData;
using ReferenceData  = S4FRTMITC::ReferenceData;

using math::normalized;

namespace {

ScalarDerivatives dot_value_B(const VectorDerivatives& a,
                              const VectorDerivatives& b) {
    ScalarDerivatives out;
    out.value = a.value.dot(b.value);
    out.d1    = a.d1 * b.value + b.d1 * a.value;
    return out;
}

void add_xx_derivatives(
    const std::array<VectorDerivatives, S4FRTMITC::num_nodes>& x_nodes,
    const StaticVector<S4FRTMITC::num_nodes>&                  a_coeffs,
    const StaticVector<S4FRTMITC::num_nodes>&                  b_coeffs,
    Precision                                                  scale,
    Index                                                      strain_id,
    Vec8&                                                      strain,
    Mat8x24*                                                   B,
    Vec24Mat*                                                  G
) {
    Vec3 a_value = Vec3::Zero();
    Vec3 b_value = Vec3::Zero();
    for (Index node = 0; node < S4FRTMITC::num_nodes; ++node) {
        a_value += a_coeffs(node) * x_nodes[node].value;
        b_value += b_coeffs(node) * x_nodes[node].value;
    }

    strain(strain_id) += scale * a_value.dot(b_value);

    if (B) {
        for (Index node = 0; node < S4FRTMITC::num_nodes; ++node) {
            const Index base = S4FRTMITC::dofs_per_node * node;
            const Vec3 value = scale
                             * (a_coeffs(node) * b_value
                              + b_coeffs(node) * a_value);

            B->template block<1, 3>(strain_id, base) += value.transpose();
        }
    }

    if (G) {
        Mat24& hessian = (*G)[strain_id];
        for (Index i = 0; i < S4FRTMITC::num_nodes; ++i) {
            const Index row = S4FRTMITC::dofs_per_node * i;
            for (Index j = 0; j < S4FRTMITC::num_nodes; ++j) {
                const Index col = S4FRTMITC::dofs_per_node * j;
                const Precision block_scale = scale
                    * (a_coeffs(i) * b_coeffs(j)
                     + b_coeffs(i) * a_coeffs(j));

                if (block_scale != Precision(0)) {
                    hessian.template block<3, 3>(row, col).diagonal().array() += block_scale;
                }
            }
        }
    }
}

void add_xd_derivatives(
    const std::array<VectorDerivatives, S4FRTMITC::num_nodes>& x_nodes,
    const std::array<VectorDerivatives, S4FRTMITC::num_nodes>& d_nodes,
    const StaticVector<S4FRTMITC::num_nodes>&                  x_coeffs,
    const StaticVector<S4FRTMITC::num_nodes>&                  d_coeffs,
    Precision                                                  scale,
    Index                                                      strain_id,
    Vec8&                                                      strain,
    Mat8x24*                                                   B,
    Vec24Mat*                                                  G
) {
    Vec3 x_value = Vec3::Zero();
    Vec3 d_value = Vec3::Zero();
    for (Index node = 0; node < S4FRTMITC::num_nodes; ++node) {
        x_value += x_coeffs(node) * x_nodes[node].value;
        d_value += d_coeffs(node) * d_nodes[node].value;
    }

    strain(strain_id) += scale * x_value.dot(d_value);

    if (B) {
        for (Index x_node = 0; x_node < S4FRTMITC::num_nodes; ++x_node) {
            const Index x_base = S4FRTMITC::dofs_per_node * x_node;
            const Vec3 value = scale * x_coeffs(x_node) * d_value;
            B->template block<1, 3>(strain_id, x_base) += value.transpose();
        }
    }

    Mat24* hessian = G ? &(*G)[strain_id] : nullptr;

    for (Index d_node = 0; d_node < S4FRTMITC::num_nodes; ++d_node) {
        const Index rot_base = S4FRTMITC::dofs_per_node * d_node + 3;
        const Precision d_coeff = d_coeffs(d_node);

        for (Index a = 0; a < 3; ++a) {
            const Vec3 director_first =
                d_nodes[d_node].d1.row(rot_base + a).transpose();

            if (B) {
                (*B)(strain_id, rot_base + a) +=
                    scale * d_coeff * director_first.dot(x_value);
            }

            if (hessian) {
                for (Index x_node = 0; x_node < S4FRTMITC::num_nodes; ++x_node) {
                    const Index x_base = S4FRTMITC::dofs_per_node * x_node;
                    const Vec3 mixed =
                        scale * x_coeffs(x_node) * d_coeff * director_first;

                    hessian->template block<3, 1>(x_base, rot_base + a) += mixed;
                    hessian->template block<1, 3>(rot_base + a, x_base) += mixed.transpose();
                }
            }
        }

        if (hessian) {
            for (Index a = 0; a < 3; ++a) {
                for (Index b = 0; b < 3; ++b) {
                    Precision second = Precision(0);
                    for (Index c = 0; c < 3; ++c) {
                        second += x_value(c)
                                * d_nodes[d_node].d2[c](rot_base + a, rot_base + b);
                    }
                    (*hessian)(rot_base + a, rot_base + b) += scale * d_coeff * second;
                }
            }
        }
    }
}

void add_xd_shear_derivatives(
    const std::array<VectorDerivatives, S4FRTMITC::num_nodes>& x_nodes,
    const std::array<VectorDerivatives, S4FRTMITC::num_nodes>& d_nodes,
    const StaticVector<S4FRTMITC::num_nodes>&                  x_coeffs,
    const StaticVector<S4FRTMITC::num_nodes>&                  d_coeffs,
    Index                                                      shear_id,
    Vec2&                                                      shear_nat,
    StaticMatrix<2, S4FRTMITC::num_dofs>*                      B_nat,
    std::array<Mat24, 2>*                                      G_nat
) {
    Vec3 x_value = Vec3::Zero();
    Vec3 d_value = Vec3::Zero();
    for (Index node = 0; node < S4FRTMITC::num_nodes; ++node) {
        x_value += x_coeffs(node) * x_nodes[node].value;
        d_value += d_coeffs(node) * d_nodes[node].value;
    }

    shear_nat(shear_id) += x_value.dot(d_value);

    if (B_nat) {
        for (Index x_node = 0; x_node < S4FRTMITC::num_nodes; ++x_node) {
            const Index x_base = S4FRTMITC::dofs_per_node * x_node;
            const Vec3 value = x_coeffs(x_node) * d_value;
            B_nat->template block<1, 3>(shear_id, x_base) += value.transpose();
        }
    }

    Mat24* hessian = G_nat ? &(*G_nat)[shear_id] : nullptr;

    for (Index d_node = 0; d_node < S4FRTMITC::num_nodes; ++d_node) {
        const Index rot_base = S4FRTMITC::dofs_per_node * d_node + 3;
        const Precision d_coeff = d_coeffs(d_node);

        for (Index a = 0; a < 3; ++a) {
            const Vec3 director_first =
                d_nodes[d_node].d1.row(rot_base + a).transpose();

            if (B_nat) {
                (*B_nat)(shear_id, rot_base + a) +=
                    d_coeff * director_first.dot(x_value);
            }

            if (hessian) {
                for (Index x_node = 0; x_node < S4FRTMITC::num_nodes; ++x_node) {
                    const Index x_base = S4FRTMITC::dofs_per_node * x_node;
                    const Vec3 mixed = x_coeffs(x_node) * d_coeff * director_first;

                    hessian->template block<3, 1>(x_base, rot_base + a) += mixed;
                    hessian->template block<1, 3>(rot_base + a, x_base) += mixed.transpose();
                }
            }
        }

        if (hessian) {
            for (Index a = 0; a < 3; ++a) {
                for (Index b = 0; b < 3; ++b) {
                    Precision second = Precision(0);
                    for (Index c = 0; c < 3; ++c) {
                        second += x_value(c)
                                * d_nodes[d_node].d2[c](rot_base + a, rot_base + b);
                    }
                    (*hessian)(rot_base + a, rot_base + b) += d_coeff * second;
                }
            }
        }
    }
}


} // namespace

S4FRTMITC::S4FRTMITC(ID id, std::array<ID, 4> nodes)
    : ShellElement<4>(id, nodes)
    , geometry(nodes)
    , integration_scheme_(quadrature::Domain::DOMAIN_ISO_QUAD,
                          quadrature::Order::ORDER_CUBIC) {}


void S4FRTMITC::step_begin() {
    auto data = std::make_unique<ReferenceData>();

    data->X    = init_ref_node_coords();
    data->area = Precision(0);

    init_reference_basis(*data);

    reference_data_ = std::move(data);
    auto& ref       = *reference_data_;

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const auto point = integration_scheme_.get_point(ip);

        init_ref_point_data(ref, ip, point.r, point.s, point.w);

        ref.area += ref.w(ip) * ref.detJ(ip);
    }

    init_ref_point_data(ref, ReferenceData::tying_start + 0, Precision(0),  Precision(-1), Precision(0));
    init_ref_point_data(ref, ReferenceData::tying_start + 1, Precision(0),  Precision( 1), Precision(0));
    init_ref_point_data(ref, ReferenceData::tying_start + 2, Precision(-1), Precision(0),  Precision(0));
    init_ref_point_data(ref, ReferenceData::tying_start + 3, Precision( 1), Precision(0),  Precision(0));
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

void S4FRTMITC::init_reference_basis(ReferenceData& ref) const {
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
    for (Index i = 0; i < num_nodes; ++i) {
        X_xi  += dN_dxi(i)  * ref.X.row(i).transpose();
        X_eta += dN_deta(i) * ref.X.row(i).transpose();
    }

    ref.e1 = normalized(X_xi);
    ref.e3 = normalized(X_xi.cross(X_eta));
    ref.e2 = normalized(ref.e3.cross(ref.e1));

    for (Index i = 0; i < num_nodes; ++i) {
        ref.d0.row(i) = ref.e3.transpose();
    }

    Vec3 x0 = Vec3::Zero();
    for (Index i = 0; i < num_nodes; ++i) {
        x0 += ref.X.row(i).transpose();
    }
    x0 /= Precision(num_nodes);

    for (Index i = 0; i < num_nodes; ++i) {
        const Vec3 dx = ref.X.row(i).transpose() - x0;
        ref.xy(i, 0) = dx.dot(ref.e1);
        ref.xy(i, 1) = dx.dot(ref.e2);
    }
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


int S4FRTMITC::cached_reference_point_id(Precision r, Precision s) const {
    if (!reference_data_) {
        return -1;
    }

    const Precision tol = Precision(1e-14);
    const auto&     ref = *reference_data_;

    for (Index ref_id = 0; ref_id < ReferenceData::num_ref_points; ++ref_id) {
        if (std::abs(ref.r(ref_id) - r) <= tol && std::abs(ref.s(ref_id) - s) <= tol) {
            return static_cast<int>(ref_id);
        }
    }

    return -1;
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

ScalarDerivatives S4FRTMITC::dot_derivatives(const VectorDerivatives& a,
                                         const VectorDerivatives& b) {
    ScalarDerivatives out;
    out.value = a.value.dot(b.value);
    out.d1    = a.d1 * b.value + b.d1 * a.value;
    out.d2.noalias() = a.d1 * b.d1.transpose();
    out.d2.noalias() += b.d1 * a.d1.transpose();

    for (Index c = 0; c < 3; ++c) {
        out.d2.noalias() += b.value(c) * a.d2[c];
        out.d2.noalias() += a.value(c) * b.d2[c];
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

void S4FRTMITC::init_ref_point_data(
    ReferenceData& ref,
    Index          ref_id,
    Precision      r,
    Precision      s,
    Precision      w
) const {
    ref.r(ref_id)    = r;
    ref.s(ref_id)    = s;
    ref.w(ref_id)    = w;
    ref.detJ(ref_id) = Precision(0);

    ref.N[ref_id]     = shape_function(r, s);
    ref.dN_rs[ref_id] = shape_derivative(r, s);
    ref.dN_da[ref_id].setZero();
    ref.dN_db[ref_id].setZero();
    ref.A[ref_id].setZero();
    ref.invA[ref_id].setZero();

    ref.X_a[ref_id].setZero();
    ref.X_b[ref_id].setZero();
    ref.X_xi[ref_id].setZero();
    ref.X_eta[ref_id].setZero();
    ref.D[ref_id].setZero();
    ref.D_a[ref_id].setZero();
    ref.D_b[ref_id].setZero();

    const StaticVector<4> a_nodes = ref.xy.col(0);
    const StaticVector<4> b_nodes = ref.xy.col(1);

    const Precision a_xi  = ref.dN_rs[ref_id].col(0).dot(a_nodes);
    const Precision a_eta = ref.dN_rs[ref_id].col(1).dot(a_nodes);
    const Precision b_xi  = ref.dN_rs[ref_id].col(0).dot(b_nodes);
    const Precision b_eta = ref.dN_rs[ref_id].col(1).dot(b_nodes);

    ref.A[ref_id] << a_xi,  b_xi,
                     a_eta, b_eta;

    const Precision detA = ref.A[ref_id].determinant();
    logging::error(
        std::abs(detA) > Precision(1e-14),
        "S4FRTMITC: singular reference Jacobian in element ",
        this->elem_id
    );

    ref.detJ(ref_id) = std::abs(detA);
    ref.invA[ref_id] = ref.A[ref_id].inverse();

    for (Index i = 0; i < num_nodes; ++i) {
        Vec2 rhs;
        rhs << ref.dN_rs[ref_id](i, 0), ref.dN_rs[ref_id](i, 1);
        const Vec2 grad = ref.invA[ref_id] * rhs;

        ref.dN_da[ref_id](i) = grad(0);
        ref.dN_db[ref_id](i) = grad(1);
    }

    for (Index i = 0; i < num_nodes; ++i) {
        const Vec3 X_i = ref.X.row(i).transpose();
        const Vec3 D_i = ref.d0.row(i).transpose();

        ref.X_a[ref_id]   += ref.dN_da[ref_id](i)    * X_i;
        ref.X_b[ref_id]   += ref.dN_db[ref_id](i)    * X_i;
        ref.X_xi[ref_id]  += ref.dN_rs[ref_id](i, 0) * X_i;
        ref.X_eta[ref_id] += ref.dN_rs[ref_id](i, 1) * X_i;
        ref.D[ref_id]     += ref.N[ref_id](i)        * D_i;
        ref.D_a[ref_id]   += ref.dN_da[ref_id](i)    * D_i;
        ref.D_b[ref_id]   += ref.dN_db[ref_id](i)    * D_i;
    }
}


CurrentState S4FRTMITC::current_state() const {
    const auto& ref = reference_data();
    const StaticMatrix<4, 6> positions = node_coords_current_6();

    CurrentState state;
    for (Index i = 0; i < 4; ++i) {
        const Vec3 x               = positions.template block<1, 3>(i, 0).transpose();
        const Vec3 rotation_vector = positions.template block<1, 3>(i, 3).transpose();
        const Vec3 d0              = ref.d0.row(i).transpose();
        const Mat3 rotation        = so3::rotation_matrix(rotation_vector);

        state.x.row(i)     = x.transpose();
        state.d.row(i)     = (rotation * d0).transpose();
        state.theta.row(i) = rotation_vector.transpose();
    }
    return state;
}


CurrentState S4FRTMITC::reference_state() const {
    const auto& ref = reference_data();
    const Mat43 X   = ref.X;

    CurrentState state;
    state.x     = X;
    state.d     = ref.d0;
    state.theta.setZero();
    return state;
}


CurrentState S4FRTMITC::current_state_from_displacement(
    const Field&        displacement) const {
    logging::error(displacement.domain == FieldDomain::NODE,
                   "S4FRTMITC: displacement field must use NODE domain");
    logging::error(displacement.components >= dofs_per_node,
                   "S4FRTMITC: displacement field requires six components");

    const auto& ref = reference_data();
    const Mat43 X   = ref.X;

    CurrentState state;
    for (Index i = 0; i < num_nodes; ++i) {
        const Index node_id = static_cast<Index>(this->node_ids[i]);
        const Vec6  q       = displacement.row_vec6(node_id);
        const Vec3  x       = X.row(i).transpose() + q.template head<3>();
        const Vec3  theta   = q.template tail<3>();
        const Vec3  d0      = ref.d0.row(i).transpose();
        const Mat3  R       = so3::rotation_matrix(theta);

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
    const auto& ref = reference_data();
    Mat3 result;
    result.col(0) = ref.e1;
    result.col(1) = ref.e2;
    result.col(2) = ref.e3;
    return result;
}


void S4FRTMITC::shape_gradients_physical(const StaticMatrix<4, 2>& dN_rs,
                              StaticVector<4>&          dN_da,
                              StaticVector<4>&          dN_db,
                              Precision&                detJ,
                              Mat2&                     A) const {
    const auto& ref = reference_data();
    const StaticVector<4> a_nodes = ref.xy.col(0);
    const StaticVector<4> b_nodes = ref.xy.col(1);

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
    const CurrentState& state,
    bool                with_second_derivatives) const {
    const auto& ref = reference_data();
    std::array<VectorDerivatives, 4> directors;

    for (Index i = 0; i < num_nodes; ++i) {
        const Index base            = dofs_per_node * i;
        const Vec3  rotation_vector = state.theta.row(i).transpose();
        const Vec3  d0              = ref.d0.row(i).transpose();

        Mat3 rotation;
        std::array<Mat3, 3> first;

        if (with_second_derivatives) {
            std::array<std::array<Mat3, 3>, 3> second;

            so3::rotation_matrix_second_derivatives(
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

            continue;
        }

        so3::rotation_matrix_first_derivatives(
            rotation_vector,
            rotation,
            first
        );

        directors[i].value = rotation * d0;

        for (Index a = 0; a < 3; ++a) {
            const Index ia = base + 3 + a;
            directors[i].d1.row(ia) = (first[a] * d0).transpose();
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
    const auto& ref = reference_data();
    const auto& X   = ref.X;

    X_a.setZero();
    X_b.setZero();
    D.setZero();
    D_a.setZero();
    D_b.setZero();

    for (Index i = 0; i < num_nodes; ++i) {
        const Vec3 X_i = X.row(i).transpose();
        const Vec3 D_i = ref.d0.row(i).transpose();

        X_a += dN_da(i) * X_i;
        X_b += dN_db(i) * X_i;
        D   += N(i)      * D_i;
        D_a += dN_da(i) * D_i;
        D_b += dN_db(i) * D_i;
    }
}
S4FRTMITC::ReferenceData::ReferenceData() {
    for (auto& value : N) {
        value.setZero();
    }
    for (auto& value : dN_rs) {
        value.setZero();
    }
    for (auto& value : dN_da) {
        value.setZero();
    }
    for (auto& value : dN_db) {
        value.setZero();
    }
    for (auto& value : A) {
        value.setZero();
    }
    for (auto& value : invA) {
        value.setZero();
    }
    for (auto& value : X_a) {
        value.setZero();
    }
    for (auto& value : X_b) {
        value.setZero();
    }
    for (auto& value : X_xi) {
        value.setZero();
    }
    for (auto& value : X_eta) {
        value.setZero();
    }
    for (auto& value : D) {
        value.setZero();
    }
    for (auto& value : D_a) {
        value.setZero();
    }
    for (auto& value : D_b) {
        value.setZero();
    }
}


S4FRTMITC::EvaluationData::EvaluationData() {
    for (auto& value : tying_shear_nat) {
        value.setZero();
    }
    for (auto& value : tying_B_nat) {
        value.setZero();
    }
    for (auto& tying : tying_G_nat) {
        for (auto& value : tying) {
            value.setZero();
        }
    }

    for (auto& value : ip_strain) {
        value.setZero();
    }
    for (auto& value : ip_B) {
        value.setZero();
    }
    for (auto& ip : ip_G) {
        for (auto& value : ip) {
            value.setZero();
        }
    }
    for (auto& value : ip_resultants) {
        value.setZero();
    }
    for (auto& value : ip_weight) {
        value = Precision(0);
    }

    eas_Kaa.setZero();
    eas_Jb.setZero();
    eas_ba.setZero();
    eas_alpha.setZero();
}


EvaluationData S4FRTMITC::init_evaluation(
    const CurrentState& state,
    bool                with_strain,
    bool                with_B,
    bool                with_G,
    bool                with_resultants,
    bool                with_eas,
    const Field*        ip_stress,
    int                 ip_start_idx
) const {
    if (with_G) {
        with_B      = true;
        with_strain = true;
    }
    if (with_B) {
        with_strain = true;
    }
    if (with_eas) {
        with_strain = true;
    }
    if (with_resultants && ip_stress == nullptr) {
        with_strain = true;
    }

    EvaluationData data;
    data.with_strain     = with_strain;
    data.with_B          = with_B;
    data.with_G          = with_G;
    data.with_resultants = with_resultants;
    data.with_eas        = with_eas;
    data.state           = state;
    data.H               = resultant_stiffness();
    data.drill_k         = drill_stiffness_per_node(data.H);

    if (with_B || with_G) {
        data.x_nodes = x_derivatives(state);
        data.d_nodes = director_derivatives(state, with_G);
    } else {
        for (Index i = 0; i < num_nodes; ++i) {
            data.x_nodes[i].value = state.x.row(i).transpose();
            data.d_nodes[i].value = state.d.row(i).transpose();
        }
    }

    if (with_strain) {
        const auto& ref = reference_data();

        for (Index tying = 0; tying < ReferenceData::num_tying_points; ++tying) {
            evaluate_tying_point(data, tying, ref, ReferenceData::tying_start + tying);
        }

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            evaluate_integration_point(data, ip, ref, ip);
        }
    }

    if (with_eas) {
        compute_eas_data(data);
    }

    if (with_resultants) {
        if (ip_stress != nullptr) {
            load_ip_resultants(data, *ip_stress, ip_start_idx);
        } else {
            compute_material_resultants(data);
        }
    }

    return data;
}


void S4FRTMITC::compute_raw_strain(
    const EvaluationData& data,
    const ReferenceData&  ref,
    Index                 ref_id,
    Vec8&                 strain,
    Mat8x24*              B,
    Vec24Mat*             G
) const {
    strain.setZero();
    if (B) {
        B->setZero();
    }
    if (G) {
        for (auto& value : *G) {
            value.setZero();
        }
    }

    const auto& N     = ref.N[ref_id];
    const auto& dN_da = ref.dN_da[ref_id];
    const auto& dN_db = ref.dN_db[ref_id];

    const bool with_G = G != nullptr;
    const bool with_B = B != nullptr || with_G;

    if (!with_B) {
        Vec3 x_a = Vec3::Zero();
        Vec3 x_b = Vec3::Zero();
        Vec3 d   = Vec3::Zero();
        Vec3 d_a = Vec3::Zero();
        Vec3 d_b = Vec3::Zero();

        for (Index i = 0; i < num_nodes; ++i) {
            x_a += dN_da(i) * data.x_nodes[i].value;
            x_b += dN_db(i) * data.x_nodes[i].value;
            d   += N(i)     * data.d_nodes[i].value;
            d_a += dN_da(i) * data.d_nodes[i].value;
            d_b += dN_db(i) * data.d_nodes[i].value;
        }

        strain(0) = Precision(0.5) * x_a.dot(x_a)
                  - Precision(0.5) * ref.X_a[ref_id].dot(ref.X_a[ref_id]);
        strain(1) = Precision(0.5) * x_b.dot(x_b)
                  - Precision(0.5) * ref.X_b[ref_id].dot(ref.X_b[ref_id]);
        strain(2) = x_a.dot(x_b)
                  - ref.X_a[ref_id].dot(ref.X_b[ref_id]);
        strain(3) = x_a.dot(d_a)
                  - ref.X_a[ref_id].dot(ref.D_a[ref_id]);
        strain(4) = x_b.dot(d_b)
                  - ref.X_b[ref_id].dot(ref.D_b[ref_id]);
        strain(5) = x_a.dot(d_b) + x_b.dot(d_a)
                  - ref.X_a[ref_id].dot(ref.D_b[ref_id])
                  - ref.X_b[ref_id].dot(ref.D_a[ref_id]);
        // we dont need these since mitc overwrites them anyway
        // strain(6) = x_a.dot(d)
        //           - ref.X_a[ref_id].dot(ref.D[ref_id]);
        // strain(7) = x_b.dot(d)
        //           - ref.X_b[ref_id].dot(ref.D[ref_id]);
        return;
    }

    if (with_G) {
        add_xx_derivatives(data.x_nodes, dN_da, dN_da, Precision(0.5), 0, strain, B, G);
        add_xx_derivatives(data.x_nodes, dN_db, dN_db, Precision(0.5), 1, strain, B, G);
        add_xx_derivatives(data.x_nodes, dN_da, dN_db, Precision(1),   2, strain, B, G);
        add_xd_derivatives(data.x_nodes, data.d_nodes, dN_da, dN_da, Precision(1), 3, strain, B, G);
        add_xd_derivatives(data.x_nodes, data.d_nodes, dN_db, dN_db, Precision(1), 4, strain, B, G);
        add_xd_derivatives(data.x_nodes, data.d_nodes, dN_da, dN_db, Precision(1), 5, strain, B, G);
        add_xd_derivatives(data.x_nodes, data.d_nodes, dN_db, dN_da, Precision(1), 5, strain, B, G);
        // we dont need these since mitc overwrites them anyway
        // add_xd_derivatives(data.x_nodes, data.d_nodes, dN_da, N,     Precision(1), 6, strain, B, G);
        // add_xd_derivatives(data.x_nodes, data.d_nodes, dN_db, N,     Precision(1), 7, strain, B, G);

        strain(0) -= Precision(0.5) * ref.X_a[ref_id].dot(ref.X_a[ref_id]);
        strain(1) -= Precision(0.5) * ref.X_b[ref_id].dot(ref.X_b[ref_id]);
        strain(2) -= ref.X_a[ref_id].dot(ref.X_b[ref_id]);
        strain(3) -= ref.X_a[ref_id].dot(ref.D_a[ref_id]);
        strain(4) -= ref.X_b[ref_id].dot(ref.D_b[ref_id]);
        strain(5) -= ref.X_a[ref_id].dot(ref.D_b[ref_id])
                   + ref.X_b[ref_id].dot(ref.D_a[ref_id]);
        // we dont need these since mitc overwrites them anyway
        // strain(6) -= ref.X_a[ref_id].dot(ref.D[ref_id]);
        // strain(7) -= ref.X_b[ref_id].dot(ref.D[ref_id]);
        return;
    }

    const VectorDerivatives x_a = linear_combination(data.x_nodes, dN_da);
    const VectorDerivatives x_b = linear_combination(data.x_nodes, dN_db);
    const VectorDerivatives d   = linear_combination(data.d_nodes, N);
    const VectorDerivatives d_a = linear_combination(data.d_nodes, dN_da);
    const VectorDerivatives d_b = linear_combination(data.d_nodes, dN_db);

    std::array<ScalarDerivatives, num_strains> items;

    items[0] = scaled(dot_value_B(x_a, x_a), Precision(0.5));
    items[1] = scaled(dot_value_B(x_b, x_b), Precision(0.5));
    items[2] = dot_value_B(x_a, x_b);
    items[3] = dot_value_B(x_a, d_a);
    items[4] = dot_value_B(x_b, d_b);

    const ScalarDerivatives xa_db = dot_value_B(x_a, d_b);
    const ScalarDerivatives xb_da = dot_value_B(x_b, d_a);

    items[5].value = xa_db.value + xb_da.value;
    items[5].d1    = xa_db.d1    + xb_da.d1;

    // we dont need these since mitc overwrites them anyway
    // items[6] = dot_value_B(x_a, d);
    // items[7] = dot_value_B(x_b, d);

    strain(0) = items[0].value - Precision(0.5) * ref.X_a[ref_id].dot(ref.X_a[ref_id]);
    strain(1) = items[1].value - Precision(0.5) * ref.X_b[ref_id].dot(ref.X_b[ref_id]);
    strain(2) = items[2].value - ref.X_a[ref_id].dot(ref.X_b[ref_id]);
    strain(3) = items[3].value - ref.X_a[ref_id].dot(ref.D_a[ref_id]);
    strain(4) = items[4].value - ref.X_b[ref_id].dot(ref.D_b[ref_id]);
    strain(5) = items[5].value
              - ref.X_a[ref_id].dot(ref.D_b[ref_id])
              - ref.X_b[ref_id].dot(ref.D_a[ref_id]);
    // we dont need these since mitc overwrites them anyway
    // strain(6) = items[6].value - ref.X_a[ref_id].dot(ref.D[ref_id]);
    // strain(7) = items[7].value - ref.X_b[ref_id].dot(ref.D[ref_id]);

    if (B) {
        for (Index a = 0; a < num_strains; ++a) {
            B->row(a) = items[a].d1.transpose();
        }
    }
}


void S4FRTMITC::compute_natural_shear(
    const EvaluationData&      data,
    const ReferenceData&       ref,
    Index                      ref_id,
    Vec2&                      shear_nat,
    StaticMatrix<2, num_dofs>* B_nat,
    std::array<Mat24, 2>*      G_nat
) const {
    shear_nat.setZero();
    if (B_nat) {
        B_nat->setZero();
    }
    if (G_nat) {
        for (auto& value : *G_nat) {
            value.setZero();
        }
    }

    const StaticVector<4> dN_dxi  = ref.dN_rs[ref_id].col(0);
    const StaticVector<4> dN_deta = ref.dN_rs[ref_id].col(1);

    const bool with_G = G_nat != nullptr;
    const bool with_B = B_nat != nullptr || with_G;

    if (!with_B) {
        Vec3 x_xi  = Vec3::Zero();
        Vec3 x_eta = Vec3::Zero();
        Vec3 d     = Vec3::Zero();

        for (Index i = 0; i < num_nodes; ++i) {
            x_xi  += dN_dxi(i)        * data.x_nodes[i].value;
            x_eta += dN_deta(i)       * data.x_nodes[i].value;
            d     += ref.N[ref_id](i) * data.d_nodes[i].value;
        }

        shear_nat(0) = x_xi.dot(d)  - ref.X_xi[ref_id].dot(ref.D[ref_id]);
        shear_nat(1) = x_eta.dot(d) - ref.X_eta[ref_id].dot(ref.D[ref_id]);
        return;
    }

    if (with_G) {
        add_xd_shear_derivatives(
            data.x_nodes,
            data.d_nodes,
            dN_dxi,
            ref.N[ref_id],
            0,
            shear_nat,
            B_nat,
            G_nat
        );
        add_xd_shear_derivatives(
            data.x_nodes,
            data.d_nodes,
            dN_deta,
            ref.N[ref_id],
            1,
            shear_nat,
            B_nat,
            G_nat
        );

        shear_nat(0) -= ref.X_xi[ref_id].dot(ref.D[ref_id]);
        shear_nat(1) -= ref.X_eta[ref_id].dot(ref.D[ref_id]);
        return;
    }

    const VectorDerivatives x_xi  = linear_combination(data.x_nodes, dN_dxi);
    const VectorDerivatives x_eta = linear_combination(data.x_nodes, dN_deta);
    const VectorDerivatives d     = linear_combination(data.d_nodes, ref.N[ref_id]);

    const ScalarDerivatives gamma_xi  = dot_value_B(x_xi, d);
    const ScalarDerivatives gamma_eta = dot_value_B(x_eta, d);

    shear_nat(0) = gamma_xi.value  - ref.X_xi[ref_id].dot(ref.D[ref_id]);
    shear_nat(1) = gamma_eta.value - ref.X_eta[ref_id].dot(ref.D[ref_id]);

    if (B_nat) {
        B_nat->row(0) = gamma_xi.d1.transpose();
        B_nat->row(1) = gamma_eta.d1.transpose();
    }
    if (G_nat) {
        (*G_nat)[0] = gamma_xi.d2;
        (*G_nat)[1] = gamma_eta.d2;
    }
}


void S4FRTMITC::compute_mitc4_shear(
    const EvaluationData&      data,
    const ReferenceData&       ref,
    Index                      ref_id,
    Vec2&                      shear,
    StaticMatrix<2, num_dofs>* B_shear,
    std::array<Mat24, 2>*      G_shear
) const {
    Vec2 g_nat = Vec2::Zero();
    StaticMatrix<2, num_dofs> B_nat;
    std::array<Mat24, 2>      G_nat;

    if (B_shear) {
        B_nat.setZero();
        B_shear->setZero();
    }
    if (G_shear) {
        for (auto& value : G_nat) {
            value.setZero();
        }
        for (auto& value : *G_shear) {
            value.setZero();
        }
    }

    const Precision r = ref.r(ref_id);
    const Precision s = ref.s(ref_id);

    const Precision c_bottom = Precision(0.5) * (Precision(1) - s);
    const Precision c_top    = Precision(0.5) * (Precision(1) + s);
    const Precision c_left   = Precision(0.5) * (Precision(1) - r);
    const Precision c_right  = Precision(0.5) * (Precision(1) + r);

    g_nat(0) = c_bottom * data.tying_shear_nat[0](0)
             + c_top    * data.tying_shear_nat[1](0);
    g_nat(1) = c_left   * data.tying_shear_nat[2](1)
             + c_right  * data.tying_shear_nat[3](1);

    if (B_shear) {
        B_nat.row(0) = c_bottom * data.tying_B_nat[0].row(0)
                     + c_top    * data.tying_B_nat[1].row(0);
        B_nat.row(1) = c_left   * data.tying_B_nat[2].row(1)
                     + c_right  * data.tying_B_nat[3].row(1);
    }

    if (G_shear) {
        G_nat[0] = c_bottom * data.tying_G_nat[0][0]
                 + c_top    * data.tying_G_nat[1][0];
        G_nat[1] = c_left   * data.tying_G_nat[2][1]
                 + c_right  * data.tying_G_nat[3][1];
    }

    shear = ref.invA[ref_id] * g_nat;

    if (B_shear) {
        *B_shear = ref.invA[ref_id] * B_nat;
    }

    if (G_shear) {
        for (Index a = 0; a < 2; ++a) {
            (*G_shear)[a].setZero();
            for (Index b = 0; b < 2; ++b) {
                (*G_shear)[a] += ref.invA[ref_id](a, b) * G_nat[b];
            }
        }
    }
}


void S4FRTMITC::evaluate_tying_point(
    EvaluationData&      data,
    Index                tying_id,
    const ReferenceData& ref,
    Index                ref_id
) const {
    compute_natural_shear(
        data,
        ref,
        ref_id,
        data.tying_shear_nat[tying_id],
        data.with_B ? &data.tying_B_nat[tying_id] : nullptr,
        data.with_G ? &data.tying_G_nat[tying_id] : nullptr
    );
}


void S4FRTMITC::evaluate_integration_point(
    EvaluationData&      data,
    Index                ip,
    const ReferenceData& ref,
    Index                ref_id
) const {
    data.ip_weight[ip] = ref.w(ref_id) * ref.detJ(ref_id);

    compute_raw_strain(
        data,
        ref,
        ref_id,
        data.ip_strain[ip],
        data.with_B ? &data.ip_B[ip] : nullptr,
        data.with_G ? &data.ip_G[ip] : nullptr
    );

    Vec2 shear;
    StaticMatrix<2, num_dofs> B_shear;
    std::array<Mat24, 2>      G_shear;

    compute_mitc4_shear(
        data,
        ref,
        ref_id,
        shear,
        data.with_B ? &B_shear : nullptr,
        data.with_G ? &G_shear : nullptr
    );

    data.ip_strain[ip](6) = shear(0);
    data.ip_strain[ip](7) = shear(1);

    if (data.with_B) {
        data.ip_B[ip].row(6) = B_shear.row(0);
        data.ip_B[ip].row(7) = B_shear.row(1);
    }

    if (data.with_G) {
        data.ip_G[ip][6] = G_shear[0];
        data.ip_G[ip][7] = G_shear[1];
    }
}


void S4FRTMITC::load_ip_resultants(
    EvaluationData& data,
    const Field&    ip_stress,
    int             ip_start_idx
) const {
    logging::error(
        ip_stress.components >= num_strains,
        "S4FRTMITC: shell resultant field requires eight components "
        "[N11,N22,N12,M11,M22,M12,Q13,Q23]"
    );

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const Index row = static_cast<Index>(ip_start_idx) + ip;
        for (Index a = 0; a < num_strains; ++a) {
            data.ip_resultants[ip](a) = ip_stress(row, a);
        }
    }
}


void S4FRTMITC::compute_eas_data(EvaluationData& data) const {
    data.eas_Kaa.setZero();
    data.eas_Jb.setZero();
    data.eas_ba.setZero();
    data.eas_alpha.setZero();

    if constexpr (eas_parameters == 0) {
        return;
    } else {
        const auto& ref = reference_data();
        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto M   = eas_matrix(ref.r(ip), ref.s(ip));
            const auto MTH = M.transpose() * data.H;
            const auto fac = data.ip_weight[ip];

            data.eas_Kaa += fac * (MTH * M);

            if (data.with_strain) {
                data.eas_ba += fac * (MTH * data.ip_strain[ip]);
            }
            if (data.with_B) {
                data.eas_Jb += fac * (MTH * data.ip_B[ip]);
            }
        }

        if (data.with_strain) {
            data.eas_alpha = -data.eas_Kaa.ldlt().solve(data.eas_ba);
        }
    }
}


void S4FRTMITC::compute_material_resultants(EvaluationData& data) const {
    logging::error(data.with_strain,
                   "S4FRTMITC: material resultants require strain evaluation");

    const auto& ref = reference_data();
    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        Vec8 eps = data.ip_strain[ip];

        if constexpr (eas_parameters > 0) {
            if (data.with_eas) {
                eps += eas_matrix(ref.r(ip), ref.s(ip)) * data.eas_alpha;
            }
        }

        data.ip_resultants[ip] = data.H * eps;
    }
}


void S4FRTMITC::assemble_material_stiffness(
    const EvaluationData& data,
    Mat24&                Kmat
) const {
    logging::error(data.with_B,
                   "S4FRTMITC: material stiffness requires B evaluation");

    Kmat.setZero();

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        Kmat += data.ip_weight[ip]
              * (data.ip_B[ip].transpose() * data.H * data.ip_B[ip]);
    }

    if constexpr (eas_parameters > 0) {
        if (data.with_eas) {
            Kmat -= data.eas_Jb.transpose()
                  * data.eas_Kaa.ldlt().solve(data.eas_Jb);
        }
    }
}


void S4FRTMITC::assemble_geometric_stiffness(
    const EvaluationData& data,
    Mat24&                Kgeo
) const {
    logging::error(data.with_G,
                   "S4FRTMITC: geometric stiffness requires G evaluation");
    logging::error(data.with_resultants,
                   "S4FRTMITC: geometric stiffness requires resultants");

    Kgeo.setZero();

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        for (Index a = 0; a < num_strains; ++a) {
            Kgeo += data.ip_weight[ip]
                  * data.ip_resultants[ip](a)
                  * data.ip_G[ip][a];
        }
    }

    Kgeo = Precision(0.5) * (Kgeo + Kgeo.transpose());
}


void S4FRTMITC::assemble_internal_force(
    const EvaluationData& data,
    Vec24&                internal_force
) const {
    logging::error(data.with_B,
                   "S4FRTMITC: internal force requires B evaluation");
    logging::error(data.with_resultants,
                   "S4FRTMITC: internal force requires resultants");

    internal_force.setZero();

    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        internal_force += data.ip_weight[ip]
                        * (data.ip_B[ip].transpose() * data.ip_resultants[ip]);
    }
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


StaticVector<S4FRTMITC::eas_parameters> S4FRTMITC::compute_eas_alpha(
    const CurrentState& state,
    const Mat8&         H
) const {
    (void) H;

    const EvaluationData data = init_evaluation(
        state,
        true,
        false,
        false,
        false,
        true
    );
    return data.eas_alpha;
}


StaticVector<S4FRTMITC::eas_parameters> S4FRTMITC::compute_eas_alpha_linear(
    const Field& displacement,
    const Mat8&  H
) const {
    StaticVector<eas_parameters> alpha;
    alpha.setZero();

    if constexpr (eas_parameters == 0) {
        (void) displacement;
        (void) H;
        return alpha;
    } else {
        const CurrentState   state = reference_state();
        const Vec24          q     = element_displacement_vector(displacement);
        const EvaluationData data  = init_evaluation(
            state,
            false,
            true,
            false,
            false,
            false
        );

        StaticMatrix<eas_parameters, eas_parameters> Kaa;
        StaticVector<eas_parameters>                 ba;
        Kaa.setZero();
        ba.setZero();

        const auto& ref = reference_data();
        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto M   = eas_matrix(ref.r(ip), ref.s(ip));
            const auto MTH = M.transpose() * H;
            const Vec8 eps = data.ip_B[ip] * q;
            const auto fac = data.ip_weight[ip];

            Kaa += fac * (MTH * M);
            ba  += fac * (MTH * eps);
        }

        alpha = -Kaa.ldlt().solve(ba);
        return alpha;
    }
}


Vec8 S4FRTMITC::generalized_strain_at(
    const Field& displacement,
    Precision    r,
    Precision    s,
    bool         nonlinear
) const {
    const Mat8 H = resultant_stiffness();

    auto make_output_ref = [&]() {
        ReferenceData out = reference_data();
        init_ref_point_data(out, 0, r, s, Precision(0));
        return out;
    };

    const int cached_id = cached_reference_point_id(r, s);

    if (nonlinear) {
        const CurrentState   state = current_state_from_displacement(displacement);
        const EvaluationData data  = init_evaluation(
            state,
            true,
            false,
            false,
            false,
            true
        );

        const ReferenceData output_ref = cached_id >= 0 ? reference_data() : make_output_ref();
        const Index         ref_id     = cached_id >= 0 ? static_cast<Index>(cached_id) : Index(0);

        Vec8 eps;
        Vec2 shear;
        compute_raw_strain(data, output_ref, ref_id, eps);
        compute_mitc4_shear(data, output_ref, ref_id, shear);
        eps(6) = shear(0);
        eps(7) = shear(1);

        if constexpr (eas_parameters > 0) {
            eps += eas_matrix(r, s) * data.eas_alpha;
        }
        return eps;
    }

    const CurrentState   state = reference_state();
    const Vec24          q     = element_displacement_vector(displacement);
    const EvaluationData data  = init_evaluation(
        state,
        false,
        true,
        false,
        false,
        false
    );

    const ReferenceData output_ref = cached_id >= 0 ? reference_data() : make_output_ref();
    const Index         ref_id     = cached_id >= 0 ? static_cast<Index>(cached_id) : Index(0);

    Vec8    dummy;
    Mat8x24 B;
    Vec2    shear;
    StaticMatrix<2, num_dofs> B_shear;

    compute_raw_strain(data, output_ref, ref_id, dummy, &B);
    compute_mitc4_shear(data, output_ref, ref_id, shear, &B_shear);

    B.row(6) = B_shear.row(0);
    B.row(7) = B_shear.row(1);

    Vec8 eps = B * q;

    if constexpr (eas_parameters > 0) {
        eps += eas_matrix(r, s) * compute_eas_alpha_linear(displacement, H);
    }
    return eps;
}


Vec8 S4FRTMITC::generalized_resultant_at(
    const Field& displacement,
    Precision    r,
    Precision    s,
    bool         nonlinear,
    Vec8*        strain_out
) const {
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


Precision S4FRTMITC::drill_stiffness_per_node(const Mat8& H) const {
    const Precision area = reference_data().area;
    const Precision shear_scale = Precision(0.5)
                                  * (std::abs(H(6, 6)) + std::abs(H(7, 7)));
    return drill_scale * shear_scale * area / Precision(num_nodes);
}


void S4FRTMITC::add_drill_stiffness(Mat24&              stiffness_matrix,
                         const Mat8&         H) const {
    const auto& ref = reference_data();
    const Precision kd = drill_stiffness_per_node(H);

    for (Index i = 0; i < num_nodes; ++i) {
        const Vec3 d0    = ref.d0.row(i).transpose();
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
    const auto& ref = reference_data();
    const Precision kd = drill_stiffness_per_node(H);

    for (Index i = 0; i < num_nodes; ++i) {
        const Vec3 d0              = ref.d0.row(i).transpose();
        const Vec3 rotation_vector = state.theta.row(i).transpose();
        const Precision drill      = rotation_vector.dot(d0);

        internal_force.template segment<3>(
            dofs_per_node * i + 3
        ) += kd * drill * d0;
    }
}
MapMatrix S4FRTMITC::stiffness(Precision* buffer) {
    const CurrentState   state = current_state();
    const EvaluationData data  = init_evaluation(
        state,
        true,
        true,
        false,
        false,
        true
    );

    Mat24 Kmat;
    assemble_material_stiffness(data, Kmat);

    add_drill_stiffness(Kmat, data.H);
    Kmat = Precision(0.5) * (Kmat + Kmat.transpose());

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = Kmat;
    return mapped;
}


MapMatrix S4FRTMITC::stiffness_geom(
    Precision*   buffer,
    const Field& ip_stress,
    int          ip_start_idx
) {
    const CurrentState   state = current_state();
    const EvaluationData data  = init_evaluation(
        state,
        false,
        false,
        true,
        true,
        false,
        &ip_stress,
        ip_start_idx
    );

    Mat24 geometric_stiffness;
    assemble_geometric_stiffness(data, geometric_stiffness);

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = geometric_stiffness;
    return mapped;
}


MapMatrix S4FRTMITC::stiffness_tangent(
    Precision*   buffer,
    NodeData&    nodal_forces,
    const Field& displacement
) {
    const CurrentState state = current_state_from_displacement(displacement);

    const EvaluationData data  = init_evaluation(
        state,
        true,
        true,
        true,
        true,
        true
    );

    Mat24 material_stiffness;
    Mat24 geometric_stiffness;
    Vec24 internal_force;

    assemble_material_stiffness  (data, material_stiffness);
    assemble_geometric_stiffness (data, geometric_stiffness);
    assemble_internal_force      (data, internal_force);

    Mat24 tangent = material_stiffness + geometric_stiffness;

    add_drill_stiffness(tangent, data.H);
    add_drill_force(internal_force, state, data.H);

    tangent = Precision(0.5) * (tangent + tangent.transpose());

    for (Index a = 0; a < num_nodes; ++a) {
        const Index node_id = static_cast<Index>(this->node_ids[a]);

        for (Index d = 0; d < dofs_per_node; ++d) {
            nodal_forces(node_id, d) += internal_force(dofs_per_node * a + d);
        }
    }

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = tangent;
    return mapped;
}


void S4FRTMITC::compute_internal_force_nonlinear(
    Field&       node_forces,
    const Field& ip_stress,
    int          ip_offset
) {
    const CurrentState   state = current_state();
    const EvaluationData data  = init_evaluation(
        state,
        false,
        true,
        false,
        true,
        false,
        &ip_stress,
        ip_offset
    );

    Vec24 internal_force;
    assemble_internal_force(data, internal_force);

    add_drill_force(internal_force, state, data.H);

    for (Index a = 0; a < num_nodes; ++a) {
        const Index node_id = static_cast<Index>(this->node_ids[a]);

        for (Index d = 0; d < dofs_per_node; ++d) {
            node_forces(node_id, d) += internal_force(dofs_per_node * a + d);
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
void S4FRTMITC::compute_stress_state(
    Field&       stress_state,
    const Field& displacement,
    int          offset,
    bool         use_green_lagrange_nl
) {
    logging::error(stress_state.components >= num_strains,
                   "S4FRTMITC: stress state requires at least eight components "
                   "[N11,N22,N12,M11,M22,M12,Q13,Q23]");

    if (use_green_lagrange_nl) {
        const CurrentState   state = current_state_from_displacement(displacement);
        const EvaluationData data  = init_evaluation(
            state,
            true,
            false,
            false,
            true,
            true
        );

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const Index row = static_cast<Index>(offset) + ip;

            for (Index component = 0; component < stress_state.components; ++component) {
                stress_state(row, component) = component < num_strains
                                             ? data.ip_resultants[ip](component)
                                             : Precision(0);
            }
        }
        return;
    }

    const Mat8           H     = resultant_stiffness();
    const CurrentState   state = reference_state();
    const EvaluationData data  = init_evaluation(
        state,
        false,
        true,
        false,
        false,
        false
    );
    const Vec24 q = element_displacement_vector(displacement);

    StaticVector<eas_parameters> alpha;
    alpha.setZero();
    if constexpr (eas_parameters > 0) {
        alpha = compute_eas_alpha_linear(displacement, H);
    }

    const auto& ref = reference_data();
    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        Vec8 eps = data.ip_B[ip] * q;

        if constexpr (eas_parameters > 0) {
            eps += eas_matrix(ref.r(ip), ref.s(ip)) * alpha;
        }

        const Vec8  values = H * eps;
        const Index row    = static_cast<Index>(offset) + ip;

        for (Index component = 0; component < stress_state.components; ++component) {
            stress_state(row, component) = component < num_strains
                                         ? values(component)
                                         : Precision(0);
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
    return this->get_section()->thickness_ * reference_data().area;
}


StaticMatrix<4, 4> S4FRTMITC::integrate_NNt() const {
    StaticMatrix<4, 4> M;
    M.setZero();

    const auto& ref = reference_data();
    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        M += ref.w(ip) * ref.detJ(ip) * (ref.N[ip] * ref.N[ip].transpose());
    }
    return M;
}


MapMatrix S4FRTMITC::mass(Precision* buffer) {
    const auto& ref = reference_data();
    const Precision    rho   = this->get_density(false);

    Mat24 mass_matrix;
    mass_matrix.setZero();

    if (rho > Precision(0)) {
        const Precision h            = this->get_section()->thickness_;
        const Precision mass_surface = rho * h;
        const Precision inertia_rot  =
            rho * h * h * h / Precision(12);
        const auto MNN = integrate_NNt();

        const Vec3 d0 = ref.e3;
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

    const auto& ref = reference_data();
    Precision result = Precision(0);
    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        result += field(global_point_current(ref.r(ip), ref.s(ip)))
                * (rho * h * ref.w(ip) * ref.detJ(ip));
    }
    return result;
}


Vec3 S4FRTMITC::integrate_vector_field(bool scale_by_density, const VecField& field) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

    const auto& ref = reference_data();
    Vec3 result = Vec3::Zero();
    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        result += field(global_point_current(ref.r(ip), ref.s(ip)))
                * (rho * h * ref.w(ip) * ref.detJ(ip));
    }
    return result;
}


void S4FRTMITC::integrate_vector_field(Field& node_loads, bool scale_by_density, const VecField& field) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

    const auto& ref = reference_data();
    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        const Vec3 f = field(global_point_current(ref.r(ip), ref.s(ip)))
                     * (rho * h * ref.w(ip) * ref.detJ(ip));

        for (Index i = 0; i < num_nodes; ++i) {
            const Index node_id = static_cast<Index>(this->node_ids[i]);
            for (Index d = 0; d < 3; ++d) {
                node_loads(node_id, d) += ref.N[ip](i) * f(d);
            }
        }
    }
}


Mat3 S4FRTMITC::integrate_tensor_field(bool scale_by_density, const TenField& field) {
    const Precision h   = this->get_section()->thickness_;
    const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

    const auto& ref = reference_data();
    Mat3 result = Mat3::Zero();
    for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
        result += field(global_point_current(ref.r(ip), ref.s(ip)))
                * (rho * h * ref.w(ip) * ref.detJ(ip));
    }
    return result;
}

void S4FRTMITC::compute_compliance(Field& displacement, Field& result) {
    Precision buffer[num_dofs * num_dofs];
    MapMatrix   K = stiffness(buffer);
    const Vec24 u = element_displacement_vector(displacement);

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
