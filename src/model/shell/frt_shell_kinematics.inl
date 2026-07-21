/**
 * @file frt_shell_kinematics.cpp
 * @brief Implements finite-rotation shell kinematics and exact derivatives.
 *
 * The implementation evaluates the current nodal directors through the SO(3)
 * exponential map and constructs the compatible Total-Lagrangian membrane,
 * bending and transverse-shear strains in natural coordinates. Their first and
 * second derivatives with respect to all element degrees of freedom are
 * evaluated analytically.
 *
 * Concrete MITC elements replace selected natural strain components after the
 * compatible evaluation. The final generalized strains and derivatives are
 * transformed into the pointwise orthonormal reference basis before the shell
 * section is evaluated.
 *
 * @see FRTShell
 *
 * @author Finn Eggers
 * @date 20.07.2026
 */

#include "frt_shell.h"

#include "../../core/logging.h"
#include "../../math/so3.h"

#include <cmath>

namespace fem::model {

namespace {

template<class VectorDerivatives, class ScalarDerivatives>
ScalarDerivatives dot_value_B(const VectorDerivatives& a,
                              const VectorDerivatives& b) {
    ScalarDerivatives result;
    result.value = a.value.dot(b.value);
    result.d1    = a.d1 * b.value + b.d1 * a.value;
    return result;
}

template<Index N>
void add_xx_derivatives(
    const std::array<typename FRTShell<N>::VectorDerivatives, N>& x_nodes,
    const typename FRTShell<N>::VecN&                            a_coefficients,
    const typename FRTShell<N>::VecN&                            b_coefficients,
    Precision                                                    scale,
    Index                                                        strain_id,
    typename FRTShell<N>::Vec8&                                  strain,
    typename FRTShell<N>::Mat8x6N*                               B,
    typename FRTShell<N>::Vec6NMat*                              G
) {
    using Shell = FRTShell<N>;

    Vec3 a_value = Vec3::Zero();
    Vec3 b_value = Vec3::Zero();

    for (Index node = 0; node < N; ++node) {
        a_value += a_coefficients(node) * x_nodes[node].value;
        b_value += b_coefficients(node) * x_nodes[node].value;
    }

    strain(strain_id) += scale * a_value.dot(b_value);

    // Differentiate the scalar product with respect to translational nodal DOFs
    if (B) {
        for (Index node = 0; node < N; ++node) {
            const Index base = Shell::dofs_per_node * node;
            const Vec3  value = scale
                              * (a_coefficients(node) * b_value
                               + b_coefficients(node) * a_value);

            B->template block<1, 3>(strain_id, base) += value.transpose();
        }
    }

    // The translational Hessian consists only of constant identity blocks
    if (G) {
        typename Shell::Mat6N& hessian = (*G)[strain_id];

        for (Index i = 0; i < N; ++i) {
            const Index row = Shell::dofs_per_node * i;

            for (Index j = 0; j < N; ++j) {
                const Index col = Shell::dofs_per_node * j;
                const Precision block_scale = scale
                    * (a_coefficients(i) * b_coefficients(j)
                     + b_coefficients(i) * a_coefficients(j));

                if (block_scale != Precision(0)) {
                    hessian.template block<3, 3>(row, col).diagonal().array() += block_scale;
                }
            }
        }
    }
}

template<Index N>
void add_xd_derivatives(
    const std::array<typename FRTShell<N>::VectorDerivatives, N>& x_nodes,
    const std::array<typename FRTShell<N>::VectorDerivatives, N>& d_nodes,
    const typename FRTShell<N>::VecN&                            x_coefficients,
    const typename FRTShell<N>::VecN&                            d_coefficients,
    Precision                                                    scale,
    Index                                                        strain_id,
    typename FRTShell<N>::Vec8&                                  strain,
    typename FRTShell<N>::Mat8x6N*                               B,
    typename FRTShell<N>::Vec6NMat*                              G
) {
    using Shell = FRTShell<N>;

    Vec3 x_value = Vec3::Zero();
    Vec3 d_value = Vec3::Zero();

    for (Index node = 0; node < N; ++node) {
        x_value += x_coefficients(node) * x_nodes[node].value;
        d_value += d_coefficients(node) * d_nodes[node].value;
    }

    strain(strain_id) += scale * x_value.dot(d_value);

    // Differentiate with respect to all translational coordinates
    if (B) {
        for (Index node = 0; node < N; ++node) {
            const Index base  = Shell::dofs_per_node * node;
            const Vec3  value = scale * x_coefficients(node) * d_value;
            B->template block<1, 3>(strain_id, base) += value.transpose();
        }
    }

    typename Shell::Mat6N* hessian = G ? &(*G)[strain_id] : nullptr;

    // Differentiate the rotated directors with respect to nodal rotations
    for (Index d_node = 0; d_node < N; ++d_node) {
        const Index     rot_base     = Shell::dofs_per_node * d_node + 3;
        const Precision d_coefficient = d_coefficients(d_node);

        for (Index a = 0; a < 3; ++a) {
            const Vec3 director_first = d_nodes[d_node].d1.row(rot_base + a).transpose();

            if (B) {
                (*B)(strain_id, rot_base + a) +=
                    scale * d_coefficient * director_first.dot(x_value);
            }

            if (hessian) {
                // Mixed translation-rotation derivative blocks
                for (Index x_node = 0; x_node < N; ++x_node) {
                    const Index x_base = Shell::dofs_per_node * x_node;
                    const Vec3  mixed  = scale
                                      * x_coefficients(x_node)
                                      * d_coefficient
                                      * director_first;

                    hessian->template block<3, 1>(x_base, rot_base + a) += mixed;
                    hessian->template block<1, 3>(rot_base + a, x_base) += mixed.transpose();
                }
            }
        }

        if (hessian) {
            // Pure rotation-rotation blocks use the second SO(3) derivatives
            for (Index a = 0; a < 3; ++a) {
                for (Index b = 0; b < 3; ++b) {
                    Precision second = Precision(0);

                    for (Index c = 0; c < 3; ++c) {
                        second += x_value(c)
                                * d_nodes[d_node].d2[c](rot_base + a, rot_base + b);
                    }

                    (*hessian)(rot_base + a, rot_base + b) +=
                        scale * d_coefficient * second;
                }
            }
        }
    }
}

template<Index N>
void transform_scalar_rows(
    const StaticMatrix<3, 3>&          transformation,
    Index                              source_start,
    Index                              target_start,
    typename FRTShell<N>::Vec8&        strain,
    typename FRTShell<N>::Mat8x6N*     B,
    typename FRTShell<N>::Vec6NMat*    G
) {
    using Shell = FRTShell<N>;

    const StaticVector<3> values = strain.template segment<3>(source_start);
    strain.template segment<3>(target_start) = transformation * values;

    if (B) {
        const StaticMatrix<3, Shell::num_dofs> rows =
            B->template block<3, Shell::num_dofs>(source_start, 0);
        B->template block<3, Shell::num_dofs>(target_start, 0) = transformation * rows;
    }

    if (G) {
        std::array<typename Shell::Mat6N, 3> source;

        for (Index component = 0; component < 3; ++component) {
            source[component] = (*G)[source_start + component];
        }

        for (Index target = 0; target < 3; ++target) {
            (*G)[target_start + target].setZero();

            for (Index source_id = 0; source_id < 3; ++source_id) {
                (*G)[target_start + target] +=
                    transformation(target, source_id) * source[source_id];
            }
        }
    }
}

} // namespace

/**
 * Collects the current global nodal positions and total rotation vectors.
 *
 * The model POSITION field stores six values per shell node in the order
 * `[x, y, z, rx, ry, rz]` for the current nonlinear configuration.
 */
template<Index N>
typename FRTShell<N>::MatN6 FRTShell<N>::node_coords_current_6() const {
    logging::error(this->_model_data            != nullptr,
                   "FRTShell: no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions != nullptr,
                   "FRTShell: POSITION field is not set");

    const auto& positions = *this->_model_data->positions;
    MatN6      values;

    for (Index node = 0; node < num_nodes; ++node) {
        values.row(node) = positions.row_vec6(static_cast<Index>(this->node_ids[node])).transpose();
    }

    return values;
}

template<Index N>
typename FRTShell<N>::CurrentState FRTShell<N>::current_state() const {
    const auto& ref       = reference_data();
    const MatN6 positions = node_coords_current_6();

    CurrentState state;

    for (Index node = 0; node < num_nodes; ++node) {
        const Vec3 x     = positions.template block<1, 3>(node, 0).transpose();
        const Vec3 theta = positions.template block<1, 3>(node, 3).transpose();
        const Vec3 d0    = ref.d0.row(node).transpose();
        const Mat3 R     = math::so3::rotation_matrix(theta);

        state.x.row(node)     = x.transpose();
        state.d.row(node)     = (R * d0).transpose();
        state.theta.row(node) = theta.transpose();
    }

    return state;
}

template<Index N>
typename FRTShell<N>::CurrentState FRTShell<N>::reference_state() const {
    const auto& ref = reference_data();

    CurrentState state;
    state.x     = ref.X;
    state.d     = ref.d0;
    state.theta.setZero();
    return state;
}

template<Index N>
typename FRTShell<N>::CurrentState FRTShell<N>::current_state_from_displacement(
    const Field& displacement
) const {
    logging::error(displacement.domain == FieldDomain::NODE,
                   "FRTShell: displacement field must use NODE domain");
    logging::error(displacement.components >= dofs_per_node,
                   "FRTShell: displacement field requires six components");

    const auto& ref = reference_data();
    CurrentState state;

    for (Index node = 0; node < num_nodes; ++node) {
        const Index node_id = static_cast<Index>(this->node_ids[node]);
        const Vec6  q       = displacement.row_vec6(node_id);
        const Vec3  x       = ref.X.row(node).transpose() + q.template head<3>();
        const Vec3  theta   = q.template tail<3>();
        const Vec3  d0      = ref.d0.row(node).transpose();
        const Mat3  R       = math::so3::rotation_matrix(theta);

        state.x.row(node)     = x.transpose();
        state.d.row(node)     = (R * d0).transpose();
        state.theta.row(node) = theta.transpose();
    }

    return state;
}

template<Index N>
typename FRTShell<N>::Vec6N FRTShell<N>::element_displacement_vector(
    const Field& displacement
) const {
    logging::error(displacement.domain == FieldDomain::NODE,
                   "FRTShell: displacement field must use NODE domain");
    logging::error(displacement.components >= dofs_per_node,
                   "FRTShell: displacement field requires six components");

    Vec6N q;

    for (Index node = 0; node < num_nodes; ++node) {
        const Index node_id = static_cast<Index>(this->node_ids[node]);
        q.template segment<6>(dofs_per_node * node) = displacement.row_vec6(node_id);
    }

    return q;
}

template<Index N>
std::array<typename FRTShell<N>::VectorDerivatives, N>
FRTShell<N>::x_derivatives(const CurrentState& state) const {
    std::array<VectorDerivatives, N> x_nodes;

    for (Index node = 0; node < num_nodes; ++node) {
        const Index base = dofs_per_node * node;

        x_nodes[node].value = state.x.row(node).transpose();
        x_nodes[node].d1(base + 0, 0) = Precision(1);
        x_nodes[node].d1(base + 1, 1) = Precision(1);
        x_nodes[node].d1(base + 2, 2) = Precision(1);
    }

    return x_nodes;
}

/**
 * Evaluates the rotated nodal directors and their exact SO(3) derivatives.
 *
 * The current director at node `i` is `d_i = R(theta_i) d0_i`. First
 * derivatives are always evaluated. Second derivatives are only requested for
 * the geometric tangent because they are considerably more expensive.
 */
template<Index N>
std::array<typename FRTShell<N>::VectorDerivatives, N>
FRTShell<N>::director_derivatives(const CurrentState& state,
                                  bool                with_second_derivatives) const {
    const auto& ref = reference_data();
    std::array<VectorDerivatives, N> directors;

    for (Index node = 0; node < num_nodes; ++node) {
        const Index base  = dofs_per_node * node;
        const Vec3  theta = state.theta.row(node).transpose();
        const Vec3  d0    = ref.d0.row(node).transpose();

        Mat3                R;
        std::array<Mat3, 3> first;

        if (with_second_derivatives) {
            std::array<std::array<Mat3, 3>, 3> second;
            math::so3::rotation_matrix_second_derivatives(theta, R, first, second);

            directors[node].value = R * d0;

            for (Index a = 0; a < 3; ++a) {
                const Index ia = base + 3 + a;
                directors[node].d1.row(ia) = (first[a] * d0).transpose();

                for (Index b = 0; b < 3; ++b) {
                    const Index ib    = base + 3 + b;
                    const Vec3  value = second[a][b] * d0;

                    for (Index component = 0; component < 3; ++component) {
                        directors[node].d2[component](ia, ib) = value(component);
                    }
                }
            }

            continue;
        }

        math::so3::rotation_matrix_first_derivatives(theta, R, first);
        directors[node].value = R * d0;

        for (Index a = 0; a < 3; ++a) {
            const Index ia = base + 3 + a;
            directors[node].d1.row(ia) = (first[a] * d0).transpose();
        }
    }

    return directors;
}

template<Index N>
typename FRTShell<N>::ScalarDerivatives FRTShell<N>::dot_derivatives(
    const VectorDerivatives& a,
    const VectorDerivatives& b
) {
    ScalarDerivatives result;
    result.value = a.value.dot(b.value);
    result.d1    = a.d1 * b.value + b.d1 * a.value;
    result.d2.noalias() = a.d1 * b.d1.transpose();
    result.d2.noalias() += b.d1 * a.d1.transpose();

    for (Index component = 0; component < 3; ++component) {
        result.d2.noalias() += b.value(component) * a.d2[component];
        result.d2.noalias() += a.value(component) * b.d2[component];
    }

    return result;
}

template<Index N>
typename FRTShell<N>::ScalarDerivatives FRTShell<N>::scaled(
    const ScalarDerivatives& value,
    Precision                scale
) {
    ScalarDerivatives result;
    result.value = scale * value.value;
    result.d1    = scale * value.d1;
    result.d2    = scale * value.d2;
    return result;
}

template<Index N>
typename FRTShell<N>::VectorDerivatives FRTShell<N>::linear_combination(
    const std::array<VectorDerivatives, N>& values,
    const VecN&                            coefficients
) {
    VectorDerivatives result;

    for (Index node = 0; node < num_nodes; ++node) {
        result.value += coefficients(node) * values[node].value;
        result.d1    += coefficients(node) * values[node].d1;

        for (Index component = 0; component < 3; ++component) {
            result.d2[component] += coefficients(node) * values[node].d2[component];
        }
    }

    return result;
}

/**
 * Evaluates the compatible generalized shell strains in natural coordinates.
 *
 * The membrane strains are Green-Lagrange strains of the current midsurface.
 * Bending strains compare the current position/director derivatives with their
 * reference values, and transverse-shear strains compare the current tangents
 * with the current interpolated director.
 *
 * The function optionally evaluates
 *
 *     B = d(epsilon) / d(q)
 *
 * and one Hessian for every generalized strain component. Engineering shear
 * conventions are used for the `rs`, curvature and transverse-shear entries.
 */
template<Index N>
void FRTShell<N>::compute_natural_strain(
    const EvaluationData& data,
    const ReferencePoint& point,
    Vec8&                 strain_nat,
    Mat8x6N*              B_nat,
    Vec6NMat*             G_nat
) const {
    using Component = ShellGeneralizedStrain::Component;

    constexpr Index epsilon_rr = static_cast<Index>(Component::EpsilonXX);
    constexpr Index epsilon_ss = static_cast<Index>(Component::EpsilonYY);
    constexpr Index gamma_rs   = static_cast<Index>(Component::GammaXY);
    constexpr Index kappa_rr   = static_cast<Index>(Component::KappaXX);
    constexpr Index kappa_ss   = static_cast<Index>(Component::KappaYY);
    constexpr Index kappa_rs   = static_cast<Index>(Component::KappaXY);
    constexpr Index gamma_r3   = static_cast<Index>(Component::GammaXZ);
    constexpr Index gamma_s3   = static_cast<Index>(Component::GammaYZ);

    strain_nat.setZero();

    if (B_nat) {
        B_nat->setZero();
    }

    if (G_nat) {
        for (auto& value : *G_nat) {
            value.setZero();
        }
    }

    const VecN dshape_r = point.dshape_rs.col(0);
    const VecN dshape_s = point.dshape_rs.col(1);
    const VecN shape     = point.shape;

    const bool with_G = G_nat != nullptr;
    const bool with_B = B_nat != nullptr || with_G;

    // Evaluate values directly when no derivatives are requested
    if (!with_B) {
        Vec3 x_r = Vec3::Zero();
        Vec3 x_s = Vec3::Zero();
        Vec3 d   = Vec3::Zero();
        Vec3 d_r = Vec3::Zero();
        Vec3 d_s = Vec3::Zero();

        for (Index node = 0; node < num_nodes; ++node) {
            x_r += dshape_r(node) * data.x_nodes[node].value;
            x_s += dshape_s(node) * data.x_nodes[node].value;
            d   += shape(node)    * data.d_nodes[node].value;
            d_r += dshape_r(node) * data.d_nodes[node].value;
            d_s += dshape_s(node) * data.d_nodes[node].value;
        }

        strain_nat(epsilon_rr) = Precision(0.5) * (x_r.dot(x_r) - point.X_r.dot(point.X_r));
        strain_nat(epsilon_ss) = Precision(0.5) * (x_s.dot(x_s) - point.X_s.dot(point.X_s));
        strain_nat(gamma_rs)   = x_r.dot(x_s) - point.X_r.dot(point.X_s);

        strain_nat(kappa_rr) = x_r.dot(d_r) - point.X_r.dot(point.D_r);
        strain_nat(kappa_ss) = x_s.dot(d_s) - point.X_s.dot(point.D_s);
        strain_nat(kappa_rs) = x_r.dot(d_s) + x_s.dot(d_r)
                             - point.X_r.dot(point.D_s) - point.X_s.dot(point.D_r);

        strain_nat(gamma_r3) = x_r.dot(d) - point.X_r.dot(point.D);
        strain_nat(gamma_s3) = x_s.dot(d) - point.X_s.dot(point.D);
        return;
    }

    // Use optimized direct assembly when second derivatives are requested
    if (with_G) {
        add_xx_derivatives<N>(data.x_nodes, dshape_r, dshape_r, Precision(0.5), epsilon_rr, strain_nat, B_nat, G_nat);
        add_xx_derivatives<N>(data.x_nodes, dshape_s, dshape_s, Precision(0.5), epsilon_ss, strain_nat, B_nat, G_nat);
        add_xx_derivatives<N>(data.x_nodes, dshape_r, dshape_s, Precision(1)  ,   gamma_rs, strain_nat, B_nat, G_nat);

        add_xd_derivatives<N>(data.x_nodes, data.d_nodes, dshape_r, dshape_r, Precision(1), kappa_rr, strain_nat, B_nat, G_nat);
        add_xd_derivatives<N>(data.x_nodes, data.d_nodes, dshape_s, dshape_s, Precision(1), kappa_ss, strain_nat, B_nat, G_nat);
        add_xd_derivatives<N>(data.x_nodes, data.d_nodes, dshape_r, dshape_s, Precision(1), kappa_rs, strain_nat, B_nat, G_nat);
        add_xd_derivatives<N>(data.x_nodes, data.d_nodes, dshape_s, dshape_r, Precision(1), kappa_rs, strain_nat, B_nat, G_nat);

        add_xd_derivatives<N>(data.x_nodes, data.d_nodes, dshape_r, shape, Precision(1), gamma_r3, strain_nat, B_nat, G_nat);
        add_xd_derivatives<N>(data.x_nodes, data.d_nodes, dshape_s, shape, Precision(1), gamma_s3, strain_nat, B_nat, G_nat);

        // Subtract the reference metric, curvature and shear contributions
        strain_nat(epsilon_rr) -= Precision(0.5) * point.X_r.dot(point.X_r);
        strain_nat(epsilon_ss) -= Precision(0.5) * point.X_s.dot(point.X_s);
        strain_nat(gamma_rs)   -= point.X_r.dot(point.X_s);
        strain_nat(kappa_rr)   -= point.X_r.dot(point.D_r);
        strain_nat(kappa_ss)   -= point.X_s.dot(point.D_s);
        strain_nat(kappa_rs)   -= point.X_r.dot(point.D_s) + point.X_s.dot(point.D_r);
        strain_nat(gamma_r3)   -= point.X_r.dot(point.D);
        strain_nat(gamma_s3)   -= point.X_s.dot(point.D);
        return;
    }

    // Construct value and first derivative from reusable vector combinations
    const VectorDerivatives x_r = linear_combination(data.x_nodes, dshape_r);
    const VectorDerivatives x_s = linear_combination(data.x_nodes, dshape_s);
    const VectorDerivatives d   = linear_combination(data.d_nodes, shape);
    const VectorDerivatives d_r = linear_combination(data.d_nodes, dshape_r);
    const VectorDerivatives d_s = linear_combination(data.d_nodes, dshape_s);

    std::array<ScalarDerivatives, num_strains> items;
    items[epsilon_rr] = scaled(dot_value_B<VectorDerivatives, ScalarDerivatives>(x_r, x_r), Precision(0.5));
    items[epsilon_ss] = scaled(dot_value_B<VectorDerivatives, ScalarDerivatives>(x_s, x_s), Precision(0.5));
    items[gamma_rs]   = dot_value_B<VectorDerivatives, ScalarDerivatives>(x_r, x_s);
    items[kappa_rr]   = dot_value_B<VectorDerivatives, ScalarDerivatives>(x_r, d_r);
    items[kappa_ss]   = dot_value_B<VectorDerivatives, ScalarDerivatives>(x_s, d_s);
    items[gamma_r3]   = dot_value_B<VectorDerivatives, ScalarDerivatives>(x_r, d);
    items[gamma_s3]   = dot_value_B<VectorDerivatives, ScalarDerivatives>(x_s, d);

    const ScalarDerivatives xr_ds = dot_value_B<VectorDerivatives, ScalarDerivatives>(x_r, d_s);
    const ScalarDerivatives xs_dr = dot_value_B<VectorDerivatives, ScalarDerivatives>(x_s, d_r);
    items[kappa_rs].value = xr_ds.value + xs_dr.value;
    items[kappa_rs].d1    = xr_ds.d1    + xs_dr.d1;

    strain_nat(epsilon_rr) = items[epsilon_rr].value - Precision(0.5) * point.X_r.dot(point.X_r);
    strain_nat(epsilon_ss) = items[epsilon_ss].value - Precision(0.5) * point.X_s.dot(point.X_s);
    strain_nat(gamma_rs)   = items[gamma_rs].value   - point.X_r.dot(point.X_s);
    strain_nat(kappa_rr)   = items[kappa_rr].value   - point.X_r.dot(point.D_r);
    strain_nat(kappa_ss)   = items[kappa_ss].value   - point.X_s.dot(point.D_s);
    strain_nat(kappa_rs)   = items[kappa_rs].value   - point.X_r.dot(point.D_s)
                                                     - point.X_s.dot(point.D_r);
    strain_nat(gamma_r3)   = items[gamma_r3].value   - point.X_r.dot(point.D);
    strain_nat(gamma_s3)   = items[gamma_s3].value   - point.X_s.dot(point.D);

    for (Index component = 0; component < num_strains; ++component) {
        B_nat->row(component) = items[component].d1.transpose();
    }
}

/**
 * Transforms generalized natural strain components into the pointwise local
 * orthonormal reference basis.
 *
 * The transformation depends exclusively on the reference geometry. It is
 * therefore applied identically to strain values, B rows and G matrices without
 * additional configuration derivatives.
 */
template<Index N>
void FRTShell<N>::transform_strain_to_local(
    const ReferencePoint& point,
    Vec8&                 strain,
    Mat8x6N*              B,
    Vec6NMat*             G
) const {
    const Precision t00 = point.invA(0, 0);
    const Precision t01 = point.invA(0, 1);
    const Precision t10 = point.invA(1, 0);
    const Precision t11 = point.invA(1, 1);

    StaticMatrix<3, 3> in_plane;
    in_plane << t00 * t00,              t01 * t01,              t00 * t01,
                t10 * t10,              t11 * t11,              t10 * t11,
                Precision(2) * t00*t10, Precision(2) * t01*t11, t00*t11 + t01*t10;

    // Transform membrane and bending tensor components
    transform_scalar_rows<N>(in_plane, 0, 0, strain, B, G);
    transform_scalar_rows<N>(in_plane, 3, 3, strain, B, G);

    // Transform the covariant transverse-shear vector
    const Vec2 shear_nat = strain.template segment<2>(6);
    strain.template segment<2>(6) = point.invA * shear_nat;

    if (B) {
        const Mat2x6N shear_B = B->template block<2, num_dofs>(6, 0);
        B->template block<2, num_dofs>(6, 0) = point.invA * shear_B;
    }

    if (G) {
        const Mat6N G_r = (*G)[6];
        const Mat6N G_s = (*G)[7];

        (*G)[6] = point.invA(0, 0) * G_r + point.invA(0, 1) * G_s;
        (*G)[7] = point.invA(1, 0) * G_r + point.invA(1, 1) * G_s;
    }
}

/**
 * Returns the engineering-strain transformation between two pointwise local
 * reference bases.
 *
 * The result maps `[eps11, eps22, gamma12]` represented in `from` into the
 * basis of `to`. Both bases belong to the undeformed reference surface, so the
 * map is constant with respect to the nonlinear degrees of freedom.
 */
template<Index N>
StaticMatrix<3, 3> FRTShell<N>::in_plane_transform(
    const ReferencePoint& from,
    const ReferencePoint& to
) const {
    const Precision c11 = to.e1.dot(from.e1);
    const Precision c12 = to.e1.dot(from.e2);
    const Precision c21 = to.e2.dot(from.e1);
    const Precision c22 = to.e2.dot(from.e2);

    StaticMatrix<3, 3> transformation;
    transformation << c11*c11,              c12*c12,              c11*c12,
                      c21*c21,              c22*c22,              c21*c22,
                      Precision(2)*c11*c21, Precision(2)*c12*c22, c11*c22 + c12*c21;
    return transformation;
}

/**
 * Evaluates all kinematic and constitutive quantities requested by one caller.
 *
 * Tying-point strains are formed first because every integration-point MITC
 * interpolation depends on them. The compatible natural strain is then
 * evaluated at each quadrature point, modified by the concrete MITC element,
 * transformed into the local reference basis and optionally passed to the
 * shell section.
 */
template<Index N>
typename FRTShell<N>::EvaluationData FRTShell<N>::init_evaluation(
    const CurrentState& state,
    bool                with_strain,
    bool                with_B,
    bool                with_G,
    bool                with_resultants,
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

    if (with_resultants && ip_stress == nullptr) {
        with_strain = true;
    }

    const auto& ref = reference_data();
    EvaluationData data(static_cast<Index>(ref.ip_points.size()),
                        static_cast<Index>(ref.tying_points.size()));

    data.with_strain     = with_strain;
    data.with_B          = with_B;
    data.with_G          = with_G;
    data.with_resultants = with_resultants;
    data.state           = state;
    data.H               = resultant_stiffness();
    data.drill_k         = drill_stiffness_per_node(data.H);

    // Initialize the integration-point section tangents for linear response and
    // for callers that do not request a fresh nonlinear material evaluation
    for (Index ip = 0; ip < static_cast<Index>(ref.ip_points.size()); ++ip) {
        const ReferencePoint& point = ref.ip_points[static_cast<std::size_t>(ip)];
        data.ip_tangent[static_cast<std::size_t>(ip)] = resultant_stiffness(point.r, point.s);
    }

    // Prepare nodal positions and directors at the requested derivative order
    if (with_B || with_G) {
        data.x_nodes = x_derivatives(state);
        data.d_nodes = director_derivatives(state, with_G);
    } else {
        for (Index node = 0; node < num_nodes; ++node) {
            data.x_nodes[node].value = state.x.row(node).transpose();
            data.d_nodes[node].value = state.d.row(node).transpose();
        }
    }

    if (with_strain) {
        evaluate_tying_points(data);
        evaluate_integration_points(data);
    }

    if (with_resultants) {
        if (ip_stress) {
            load_ip_resultants(data, *ip_stress, ip_start_idx);
        } else {
            compute_material_resultants(data);
        }
    }

    return data;
}

template<Index N>
void FRTShell<N>::evaluate_tying_points(EvaluationData& data) const {
    const auto& points = reference_data().tying_points;

    for (Index tying = 0; tying < static_cast<Index>(points.size()); ++tying) {
        compute_natural_strain(
            data,
            points[static_cast<std::size_t>(tying)],
            data.tying_strain_nat[static_cast<std::size_t>(tying)],
            data.with_B ? &data.tying_B_nat[static_cast<std::size_t>(tying)] : nullptr,
            data.with_G ? &data.tying_G_nat[static_cast<std::size_t>(tying)] : nullptr
        );
    }
}

template<Index N>
void FRTShell<N>::evaluate_integration_points(EvaluationData& data) const {
    const auto& points = reference_data().ip_points;

    for (Index ip = 0; ip < static_cast<Index>(points.size()); ++ip) {
        const ReferencePoint& point = points[static_cast<std::size_t>(ip)];
        Vec8&                 strain = data.ip_strain[static_cast<std::size_t>(ip)];
        Mat8x6N* B = data.with_B ? &data.ip_B[static_cast<std::size_t>(ip)] : nullptr;
        Vec6NMat* G = data.with_G ? &data.ip_G[static_cast<std::size_t>(ip)] : nullptr;

        // Evaluate the compatible natural strain field
        compute_natural_strain(data, point, strain, B, G);

        // Replace the selected components by the element-specific MITC field
        apply_mitc_natural(data, point, strain, B, G);

        // Express all generalized quantities in the local material basis
        transform_strain_to_local(point, strain, B, G);

        data.ip_weight[static_cast<std::size_t>(ip)] = point.w * point.detJ;
    }
}

} // namespace fem::model
