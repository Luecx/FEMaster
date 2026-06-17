/**
 * @file s4_frt_mitc_nonlinear.h
 * @brief Four-node fully geometrically nonlinear FRT/MITC shell element.
 *
 * The element is a Total-Lagrangian Reissner-Mindlin shell with finite nodal
 * rotations. Each solver node contributes six global DOFs ordered as
 * [ux, uy, uz, rx, ry, rz]. The nonlinear solver stores the physical total
 * orientation as a logarithmic rotation vector and updates the orientation
 * multiplicatively on SO(3). The element linearization is therefore written
 * with respect to spatial incremental rotations about the current state.
 *
 * Mechanical model:
 *   - Total-Lagrangian reference-midsurface integration.
 *   - Current director d = R d0 with R = Exp(theta).
 *   - Green-Lagrange membrane, bending and transverse-shear strains.
 *   - MITC4 transverse-shear tying, including differentiated B and G terms.
 *   - Optional compact membrane EAS condensation controlled by eas_parameters.
 *   - Split nonlinear tangent K_T = K_mat + K_geo.
 *   - K_geo contains membrane, bending and transverse-shear resultant terms.
 *
 * A small tangent-only drilling stabilization regularizes the rotation about
 * the current director. It is a numerical stabilization, not an additional
 * physical shell resultant.
 *
 * The reference director and in-plane basis are generated locally from the
 * planar reference element geometry.
 */

#pragma once

#include "shell.h"
#include "../geometry/surface/surface4.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>

namespace fem::model {

struct S4FRTMITC : ShellElement<4> {
    using Base = ShellElement<4>;

    static constexpr Index num_nodes      = 4;
    static constexpr Index dofs_per_node = 6;
    static constexpr Index num_dofs       = num_nodes * dofs_per_node;
    static constexpr Index num_strains    = 8;

    // 0 disables EAS. Use 4 for the compact membrane EAS4 mode or 5 for EAS5.
    static constexpr Index eas_parameters = 0;

    // Small tangent-only drill stabilization. It regularizes the otherwise
    // rank-deficient rotation about the current director without introducing a
    // physical drilling energy into the shell resultants.
    static constexpr Precision drill_scale = Precision(1e-5);

    using Vec8     = StaticVector<num_strains>;
    using Vec24    = StaticVector<num_dofs>;
    using Mat8     = StaticMatrix<num_strains, num_strains>;
    using Mat24    = StaticMatrix<num_dofs, num_dofs>;
    using Mat8x24  = StaticMatrix<num_strains, num_dofs>;
    using Mat43    = StaticMatrix<4, 3>;
    using Mat42    = StaticMatrix<4, 2>;
    using Vec24Mat = std::array<Mat24, num_strains>;

    struct VectorDerivatives {
        Vec3                        value = Vec3::Zero();
        StaticMatrix<num_dofs, 3>   d1    = StaticMatrix<num_dofs, 3>::Zero();
        std::array<Mat24, 3>        d2{};

        VectorDerivatives() {
            for (auto& matrix : d2) {
                matrix.setZero();
            }
        }
    };

    struct ScalarDerivatives {
        Precision value = Precision(0);
        Vec24     d1    = Vec24::Zero();
        Mat24     d2    = Mat24::Zero();
    };

    struct ElementBasis {
        Vec3  e1 = Vec3::Zero();
        Vec3  e2 = Vec3::Zero();
        Vec3  e3 = Vec3::Zero();
        Mat43 d0 = Mat43::Zero();
        Mat42 xy = Mat42::Zero();
    };

    struct CurrentState {
        Mat43                 x = Mat43::Zero();
        Mat43                 d = Mat43::Zero();
        std::array<Mat3, 4>   rotation{};
    };

    struct StrainData {
        Vec8     strain = Vec8::Zero();
        Mat8x24  B      = Mat8x24::Zero();
        Vec24Mat G{};
        Precision detJ  = Precision(0);

        StrainData() {
            for (auto& matrix : G) {
                matrix.setZero();
            }
        }
    };

    Surface4                 geometry;
    quadrature::Quadrature   integration_scheme_;

    S4FRTMITC(ID id, std::array<ID, 4> nodes)
        : ShellElement<4>(id, nodes)
        , geometry(nodes)
        , integration_scheme_(quadrature::Domain::DOMAIN_ISO_QUAD,
                              quadrature::Order::ORDER_CUBIC) {}

    ~S4FRTMITC() override = default;

    std::string type_name() const override { return "MITC4FRT"; }

    std::shared_ptr<SurfaceInterface> surface(int surface_id) override {
        return std::make_shared<Surface4>(
            surface_id == 1
                ? std::array<ID, 4>{this->nodes()[0], this->nodes()[1], this->nodes()[2], this->nodes()[3]}
                : std::array<ID, 4>{this->nodes()[3], this->nodes()[2], this->nodes()[1], this->nodes()[0]}
        );
    }

    const quadrature::Quadrature& integration_scheme() const override {
        return integration_scheme_;
    }

    StaticVector<4> shape_function(Precision r, Precision s) const {
        StaticVector<4> N;
        N(0) = Precision(0.25) * (Precision(1) - r) * (Precision(1) - s);
        N(1) = Precision(0.25) * (Precision(1) + r) * (Precision(1) - s);
        N(2) = Precision(0.25) * (Precision(1) + r) * (Precision(1) + s);
        N(3) = Precision(0.25) * (Precision(1) - r) * (Precision(1) + s);
        return N;
    }

    StaticMatrix<4, 2> shape_derivative(Precision r, Precision s) const {
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

    ShellSection* shell_section() const {
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

    static Vec3 normalized(Vec3 v, const std::string& name) {
        const Precision n = v.norm();
        logging::error(n > Precision(1e-14), "S4FRTMITC: cannot normalize ", name);
        return v / n;
    }

    static Mat3 skew(const Vec3& v) {
        Mat3 S;
        S << Precision(0), -v(2),       v(1),
             v(2),        Precision(0), -v(0),
            -v(1),        v(0),        Precision(0);
        return S;
    }

    static Mat3 rotation_exp(const Vec3& rotation_vector) {
        const Precision angle_squared = rotation_vector.squaredNorm();
        const Precision angle         = std::sqrt(angle_squared);
        const Mat3      omega         = skew(rotation_vector);
        const Mat3      omega_squared = omega * omega;

        Precision a;
        Precision b;

        if (angle < Precision(1e-8)) {
            const Precision angle_fourth = angle_squared * angle_squared;
            a = Precision(1)
              - angle_squared / Precision(6)
              + angle_fourth  / Precision(120);
            b = Precision(0.5)
              - angle_squared / Precision(24)
              + angle_fourth  / Precision(720);
        } else {
            a = std::sin(angle) / angle;
            b = (Precision(1) - std::cos(angle)) / angle_squared;
        }

        return Mat3::Identity() + a * omega + b * omega_squared;
    }

    static ScalarDerivatives dot_derivatives(const VectorDerivatives& a,
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

    static ScalarDerivatives scaled(const ScalarDerivatives& a, Precision scale) {
        ScalarDerivatives out;
        out.value = scale * a.value;
        out.d1    = scale * a.d1;
        out.d2    = scale * a.d2;
        return out;
    }

    static VectorDerivatives linear_combination(const std::array<VectorDerivatives, 4>& values,
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

    Mat43 node_coords_reference_xyz() const {
        logging::error(this->_model_data != nullptr,
                       "S4FRTMITC: no model data assigned to element ", this->elem_id);
        logging::error(this->_model_data->positions_reference != nullptr,
                       "S4FRTMITC: POSITION_REFERENCE field is not set");

        const auto& positions = *this->_model_data->positions_reference;
        Mat43 X;
        for (Index i = 0; i < 4; ++i) {
            X.row(i) = positions.row_vec3(static_cast<Index>(this->node_ids[i])).transpose();
        }
        return X;
    }

    StaticMatrix<4, 6> node_coords_current_6() const {
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

    ElementBasis reference_basis() const {
        const Mat43 X = node_coords_reference_xyz();
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

        basis.e1 = normalized(X_xi, "reference axis 1");
        basis.e3 = normalized(X_xi.cross(X_eta), "reference director");
        basis.e2 = normalized(basis.e3.cross(basis.e1), "reference axis 2");

        for (Index i = 0; i < 4; ++i) {
            basis.d0.row(i) = basis.e3.transpose();
        }

        compute_local_reference_coordinates(X, basis);
        return basis;
    }

    static void compute_local_reference_coordinates(const Mat43& X, ElementBasis& basis) {
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
    }

    CurrentState current_state(const ElementBasis& basis) const {
        const StaticMatrix<4, 6> positions = node_coords_current_6();

        CurrentState state;
        for (Index i = 0; i < 4; ++i) {
            const Vec3 x               = positions.template block<1, 3>(i, 0).transpose();
            const Vec3 rotation_vector = positions.template block<1, 3>(i, 3).transpose();
            const Vec3 d0              = basis.d0.row(i).transpose();
            const Mat3 rotation        = rotation_exp(rotation_vector);

            state.x.row(i)       = x.transpose();
            state.d.row(i)       = (rotation * d0).transpose();
            state.rotation[i]    = rotation;
        }
        return state;
    }

    void shape_gradients_physical(const StaticMatrix<4, 2>& dN_rs,
                                  const ElementBasis&       basis,
                                  StaticVector<4>&          dN_da,
                                  StaticVector<4>&          dN_db,
                                  Precision&                detJ,
                                  Mat2&                     A) const {
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

    std::array<VectorDerivatives, 4> x_derivatives(const CurrentState& state) const {
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

    std::array<VectorDerivatives, 4> director_derivatives(const CurrentState& state) const {
        std::array<VectorDerivatives, 4> directors;

        std::array<Vec3, 3> axes{
            Vec3::UnitX(),
            Vec3::UnitY(),
            Vec3::UnitZ()
        };

        std::array<Mat3, 3> generators{
            skew(axes[0]),
            skew(axes[1]),
            skew(axes[2])
        };

        for (Index i = 0; i < 4; ++i) {
            const Index base     = dofs_per_node * i;
            const Vec3  director = state.d.row(i).transpose();

            directors[i].value = director;

            // The nonlinear solver supplies spatial incremental rotations and
            // updates R_{k+1} = Exp(delta_theta) R_k. Therefore the local
            // derivatives are evaluated at zero increment about the current
            // orientation, not with respect to the components of Log(R).
            for (Index a = 0; a < 3; ++a) {
                const Index ia = base + 3 + a;
                directors[i].d1.row(ia) = (generators[a] * director).transpose();

                for (Index b = 0; b < 3; ++b) {
                    const Index ib = base + 3 + b;
                    const Vec3 second = Precision(0.5)
                                      * (generators[a] * generators[b]
                                       + generators[b] * generators[a])
                                      * director;

                    for (Index c = 0; c < 3; ++c) {
                        directors[i].d2[c](ia, ib) = second(c);
                    }
                }
            }
        }

        return directors;
    }

    void reference_fields(const ElementBasis& basis,
                          const StaticVector<4>& N,
                          const StaticVector<4>& dN_da,
                          const StaticVector<4>& dN_db,
                          Vec3& X_a,
                          Vec3& X_b,
                          Vec3& D,
                          Vec3& D_a,
                          Vec3& D_b) const {
        const Mat43 X = node_coords_reference_xyz();

        X_a.setZero();
        X_b.setZero();
        D.setZero();
        D_a.setZero();
        D_b.setZero();

        for (Index i = 0; i < 4; ++i) {
            X_a += dN_da(i) * X.row(i).transpose();
            X_b += dN_db(i) * X.row(i).transpose();
            D   += N(i)      * basis.d0.row(i).transpose();
            D_a += dN_da(i) * basis.d0.row(i).transpose();
            D_b += dN_db(i) * basis.d0.row(i).transpose();
        }
    }

    StrainData raw_strain_B_G_at(const CurrentState& state,
                                 const ElementBasis& basis,
                                 Precision           r,
                                 Precision           s) const {
        StrainData out;

        const StaticVector<4>   N     = shape_function(r, s);
        const StaticMatrix<4,2> dN_rs = shape_derivative(r, s);

        StaticVector<4> dN_da;
        StaticVector<4> dN_db;
        Mat2            A;
        shape_gradients_physical(dN_rs, basis, dN_da, dN_db, out.detJ, A);

        const auto x_nodes = x_derivatives(state);
        const auto d_nodes = director_derivatives(state);

        Vec3 X_a, X_b, D, D_a, D_b;
        reference_fields(basis, N, dN_da, dN_db, X_a, X_b, D, D_a, D_b);

        const VectorDerivatives x_a = linear_combination(x_nodes, dN_da);
        const VectorDerivatives x_b = linear_combination(x_nodes, dN_db);
        const VectorDerivatives d   = linear_combination(d_nodes, N);
        const VectorDerivatives d_a = linear_combination(d_nodes, dN_da);
        const VectorDerivatives d_b = linear_combination(d_nodes, dN_db);

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
            Precision(0.5) * X_a.dot(X_a),
            Precision(0.5) * X_b.dot(X_b),
            X_a.dot(X_b),
            X_a.dot(D_a),
            X_b.dot(D_b),
            X_a.dot(D_b) + X_b.dot(D_a),
            X_a.dot(D),
            X_b.dot(D)
        }};

        for (Index a = 0; a < 8; ++a) {
            out.strain(a) = items[a].value - constants[a];
            out.B.row(a)  = items[a].d1.transpose();
            out.G[a]      = items[a].d2;
        }

        return out;
    }

    void raw_natural_shear_B_G_at(const CurrentState& state,
                                  const ElementBasis& basis,
                                  Precision           r,
                                  Precision           s,
                                  Vec2&               shear_nat,
                                  StaticMatrix<2, num_dofs>& B_nat,
                                  std::array<Mat24, 2>& G_nat) const {
        const StaticVector<4>   N     = shape_function(r, s);
        const StaticMatrix<4,2> dN_rs = shape_derivative(r, s);

        const auto x_nodes = x_derivatives(state);
        const auto d_nodes = director_derivatives(state);

        StaticVector<4> dN_dxi  = dN_rs.col(0);
        StaticVector<4> dN_deta = dN_rs.col(1);

        const VectorDerivatives x_xi  = linear_combination(x_nodes, dN_dxi);
        const VectorDerivatives x_eta = linear_combination(x_nodes, dN_deta);
        const VectorDerivatives d     = linear_combination(d_nodes, N);

        const ScalarDerivatives gamma_xi  = dot_derivatives(x_xi, d);
        const ScalarDerivatives gamma_eta = dot_derivatives(x_eta, d);

        const Mat43 X = node_coords_reference_xyz();
        Vec3 X_xi  = Vec3::Zero();
        Vec3 X_eta = Vec3::Zero();
        Vec3 D     = Vec3::Zero();

        for (Index i = 0; i < 4; ++i) {
            X_xi  += dN_dxi(i)  * X.row(i).transpose();
            X_eta += dN_deta(i) * X.row(i).transpose();
            D     += N(i)       * basis.d0.row(i).transpose();
        }

        shear_nat(0) = gamma_xi.value  - X_xi.dot(D);
        shear_nat(1) = gamma_eta.value - X_eta.dot(D);

        B_nat.row(0) = gamma_xi.d1.transpose();
        B_nat.row(1) = gamma_eta.d1.transpose();

        G_nat[0] = gamma_xi.d2;
        G_nat[1] = gamma_eta.d2;
    }

    void mitc4_shear_B_G_at(const CurrentState& state,
                            const ElementBasis& basis,
                            Precision           r,
                            Precision           s,
                            const Mat2&         A,
                            Vec2&               shear,
                            StaticMatrix<2, num_dofs>& B_shear,
                            std::array<Mat24, 2>& G_shear) const {
        Vec2 gb, gt, gl, gr;
        StaticMatrix<2, num_dofs> Bb, Bt, Bl, Br;
        std::array<Mat24, 2> Gb, Gt, Gl, Gr;
        for (Index i = 0; i < 2; ++i) {
            Gb[i].setZero(); Gt[i].setZero(); Gl[i].setZero(); Gr[i].setZero();
        }

        raw_natural_shear_B_G_at(state, basis, Precision(0), Precision(-1), gb, Bb, Gb);
        raw_natural_shear_B_G_at(state, basis, Precision(0), Precision( 1), gt, Bt, Gt);
        raw_natural_shear_B_G_at(state, basis, Precision(-1), Precision(0), gl, Bl, Gl);
        raw_natural_shear_B_G_at(state, basis, Precision( 1), Precision(0), gr, Br, Gr);

        Vec2 g_nat = Vec2::Zero();
        StaticMatrix<2, num_dofs> B_nat = StaticMatrix<2, num_dofs>::Zero();
        std::array<Mat24, 2> G_nat;
        for (auto& m : G_nat) {
            m.setZero();
        }

        const Precision c_bottom = Precision(0.5) * (Precision(1) - s);
        const Precision c_top    = Precision(0.5) * (Precision(1) + s);
        const Precision c_left   = Precision(0.5) * (Precision(1) - r);
        const Precision c_right  = Precision(0.5) * (Precision(1) + r);

        g_nat(0)    = c_bottom * gb(0)      + c_top   * gt(0);
        B_nat.row(0)= c_bottom * Bb.row(0)  + c_top   * Bt.row(0);
        G_nat[0]    = c_bottom * Gb[0]      + c_top   * Gt[0];

        g_nat(1)    = c_left   * gl(1)      + c_right * gr(1);
        B_nat.row(1)= c_left   * Bl.row(1)  + c_right * Br.row(1);
        G_nat[1]    = c_left   * Gl[1]      + c_right * Gr[1];

        const Mat2 invA = A.inverse();
        shear           = invA * g_nat;
        B_shear         = invA * B_nat;
        for (Index a = 0; a < 2; ++a) {
            G_shear[a].setZero();
            for (Index b = 0; b < 2; ++b) {
                G_shear[a] += invA(a, b) * G_nat[b];
            }
        }
    }

    StrainData strain_B_G_at(const CurrentState& state,
                             const ElementBasis& basis,
                             Precision           r,
                             Precision           s) const {
        StrainData out = raw_strain_B_G_at(state, basis, r, s);

        const StaticMatrix<4,2> dN_rs = shape_derivative(r, s);
        StaticVector<4> dN_da;
        StaticVector<4> dN_db;
        Precision       detJ_unused;
        Mat2            A;
        shape_gradients_physical(dN_rs, basis, dN_da, dN_db, detJ_unused, A);

        Vec2 shear;
        StaticMatrix<2, num_dofs> B_shear;
        std::array<Mat24, 2> G_shear;
        mitc4_shear_B_G_at(state, basis, r, s, A, shear, B_shear, G_shear);

        out.strain(6) = shear(0);
        out.strain(7) = shear(1);
        out.B.row(6)  = B_shear.row(0);
        out.B.row(7)  = B_shear.row(1);
        out.G[6]      = G_shear[0];
        out.G[7]      = G_shear[1];

        return out;
    }

    // ---------------------------------------------------------------------
    // Material/resultant matrices and EAS
    // ---------------------------------------------------------------------

    Mat8 resultant_stiffness() const {
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

    StaticMatrix<num_strains, eas_parameters> eas_matrix(Precision r, Precision s) const {
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

    StaticVector<eas_parameters> compute_eas_alpha(const CurrentState& state,
                                                   const ElementBasis& basis,
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

            for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
                const auto      p = integration_scheme_.get_point(ip);
                const StrainData sd = strain_B_G_at(state, basis, p.r, p.s);
                const Precision fac = p.w * sd.detJ;
                const auto      M   = eas_matrix(p.r, p.s);
                const StaticMatrix<eas_parameters, num_strains> MTH = M.transpose() * H;

                Kaa += fac * (MTH * M);
                ba  += fac * (MTH * sd.strain);
            }

            alpha = -Kaa.ldlt().solve(ba);
            return alpha;
        }
    }

    void material_and_geometric_stiffness(const CurrentState& state,
                                            const ElementBasis& basis,
                                            const Mat8&         H,
                                            Mat24&              Kmat,
                                            Mat24&              Kgeo,
                                            Vec24*              force = nullptr) const {
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

            for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
                const auto       p   = integration_scheme_.get_point(ip);
                const StrainData sd  = strain_B_G_at(state, basis, p.r, p.s);
                const Precision  fac = p.w * sd.detJ;
                const Vec8       res = H * sd.strain;

                Kmat += fac * (sd.B.transpose() * H * sd.B);
                for (Index a = 0; a < num_strains; ++a) {
                    Kgeo += fac * res(a) * sd.G[a];
                }
                if (force) {
                    *force += fac * (sd.B.transpose() * res);
                }

                const auto M = eas_matrix(p.r, p.s);
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
            for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
                const auto      p  = integration_scheme_.get_point(ip);
                const StrainData sd = strain_B_G_at(state, basis, p.r, p.s);
                const Precision fac = p.w * sd.detJ;
                const Vec8      s  = H * sd.strain;

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

    Vec8 shell_resultant_at(const CurrentState& state,
                            const ElementBasis& basis,
                            const Mat8&         H,
                            Precision           r,
                            Precision           s,
                            Vec8*               strain_out = nullptr) const {
        const StrainData sd = strain_B_G_at(state, basis, r, s);
        Vec8 eps = sd.strain;

        if constexpr (eas_parameters > 0) {
            const auto alpha = compute_eas_alpha(state, basis, H);
            eps += eas_matrix(r, s) * alpha;
        }

        if (strain_out) {
            *strain_out = eps;
        }
        return H * eps;
    }

    // ---------------------------------------------------------------------
    // Drill stabilization and StructuralElement interface
    // ---------------------------------------------------------------------

    Precision reference_area(const ElementBasis& basis) const {
        Precision area = Precision(0);

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto      point = integration_scheme_.get_point(ip);
            const auto      dN_rs = shape_derivative(point.r, point.s);
            StaticVector<4> dN_da;
            StaticVector<4> dN_db;
            Precision       detJ;
            Mat2            jacobian;

            shape_gradients_physical(
                dN_rs,
                basis,
                dN_da,
                dN_db,
                detJ,
                jacobian
            );

            area += point.w * detJ;
        }

        return area;
    }

    Precision drill_stiffness_per_node(const ElementBasis& basis,
                                       const Mat8&         H) const {
        const Precision area = reference_area(basis);
        const Precision shear_scale = Precision(0.5)
                                      * (std::abs(H(6, 6)) + std::abs(H(7, 7)));
        return drill_scale * shear_scale * area / Precision(num_nodes);
    }

    void add_drill_stiffness(Mat24&              stiffness_matrix,
                             const CurrentState& state,
                             const ElementBasis& basis,
                             const Mat8&         H) const {
        const Precision kd = drill_stiffness_per_node(basis, H);

        for (Index i = 0; i < num_nodes; ++i) {
            const Vec3 director = state.d.row(i).transpose();
            const Mat3 block    = kd * (director * director.transpose());

            stiffness_matrix.template block<3, 3>(
                dofs_per_node * i + 3,
                dofs_per_node * i + 3
            ) += block;
        }
    }

    MapMatrix stiffness(Precision* buffer) override {
        const ElementBasis basis = reference_basis();
        const CurrentState state = current_state(basis);
        const Mat8         H     = resultant_stiffness();

        Mat24 material_stiffness;
        Mat24 geometric_dummy;

        material_and_geometric_stiffness(
            state,
            basis,
            H,
            material_stiffness,
            geometric_dummy,
            nullptr
        );

        add_drill_stiffness(material_stiffness, state, basis, H);
        material_stiffness = Precision(0.5)
                           * (material_stiffness + material_stiffness.transpose());

        MapMatrix mapped(buffer, num_dofs, num_dofs);
        mapped = material_stiffness;
        return mapped;
    }

    MapMatrix stiffness_geom(Precision*   buffer,
                             const Field& ip_stress,
                             int          ip_start_idx) override {
        logging::error(
            ip_stress.components >= num_strains,
            "S4FRTMITC: geometric stiffness requires eight shell resultants "
            "[N11,N22,N12,M11,M22,M12,Q13,Q23]"
        );

        const ElementBasis basis = reference_basis();
        const CurrentState state = current_state(basis);

        Mat24 geometric_stiffness;
        geometric_stiffness.setZero();

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto       point  = integration_scheme_.get_point(ip);
            const StrainData data   = strain_B_G_at(state, basis, point.r, point.s);
            const Precision  factor = point.w * data.detJ;
            const Index      row    = static_cast<Index>(ip_start_idx) + ip;

            Vec8 resultants;
            for (Index a = 0; a < num_strains; ++a) {
                resultants(a) = ip_stress(row, a);
            }

            // Rao's finite-rotation tangent contracts the complete generalized
            // resultant vector with the second strain variation. Consequently,
            // membrane, bending and transverse-shear contributions all enter.
            for (Index a = 0; a < num_strains; ++a) {
                geometric_stiffness += factor * resultants(a) * data.G[a];
            }
        }

        geometric_stiffness = Precision(0.5)
                            * (geometric_stiffness + geometric_stiffness.transpose());

        MapMatrix mapped(buffer, num_dofs, num_dofs);
        mapped = geometric_stiffness;
        return mapped;
    }

    void compute_internal_force_nonlinear(Field&       node_forces,
                                          const Field& ip_stress,
                                          int          ip_offset) override {
        logging::error(
            ip_stress.components >= num_strains,
            "S4FRTMITC: internal force requires eight shell resultants "
            "[N11,N22,N12,M11,M22,M12,Q13,Q23]"
        );

        const ElementBasis basis = reference_basis();
        const CurrentState state = current_state(basis);

        Vec24 internal_force = Vec24::Zero();

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto       point  = integration_scheme_.get_point(ip);
            const StrainData data   = strain_B_G_at(state, basis, point.r, point.s);
            const Precision  factor = point.w * data.detJ;
            const Index      row    = static_cast<Index>(ip_offset) + ip;

            Vec8 resultants;
            for (Index a = 0; a < num_strains; ++a) {
                resultants(a) = ip_stress(row, a);
            }

            internal_force += factor * (data.B.transpose() * resultants);
        }

        // The drilling stabilization is tangent-only. It regularizes the null
        // mode but does not add a physical drilling resultant.
        for (Index a = 0; a < num_nodes; ++a) {
            const Index node_id = static_cast<Index>(this->node_ids[a]);

            for (Index d = 0; d < dofs_per_node; ++d) {
                node_forces(node_id, d) +=
                    internal_force(dofs_per_node * a + d);
            }
        }
    }

    void compute_stress_strain(Field*           strain,
                               Field*           stress_out,
                               const Field&     displacement,
                               const RowMatrix& rst,
                               int              offset,
                               bool             use_green_lagrange_nl) override {
        (void) displacement;
        (void) use_green_lagrange_nl;

        logging::error(strain != nullptr || stress_out != nullptr,
                       "S4FRTMITC: compute_stress_strain requires at least one output field");
        logging::error(rst.cols() >= 2,
                       "S4FRTMITC: stress/strain coordinates require at least r,s columns");

        const ElementBasis basis = reference_basis();
        const CurrentState state = current_state(basis);
        const Mat8         H     = resultant_stiffness();

        for (Index n = 0; n < rst.rows(); ++n) {
            const Precision r   = rst(n, 0);
            const Precision s   = rst(n, 1);
            const Index     row = static_cast<Index>(offset) + n;

            Vec8 eps;
            const Vec8 res = shell_resultant_at(state, basis, H, r, s, &eps);

            if (strain) {
                const Index ncomp = std::min<Index>(strain->components, num_strains);
                for (Index a = 0; a < ncomp; ++a) {
                    (*strain)(row, a) = eps(a);
                }
            }
            if (stress_out) {
                const Index ncomp = std::min<Index>(stress_out->components, num_strains);
                for (Index a = 0; a < ncomp; ++a) {
                    (*stress_out)(row, a) = res(a);
                }
            }
        }
    }

    RowMatrix stress_strain_ip_rst() override {
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

    RowMatrix stress_strain_nodal_rst() override {
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

    Precision volume() override {
        const ElementBasis basis = reference_basis();
        return this->get_section()->thickness_ * reference_area(basis);
    }

    StaticMatrix<4, 4> integrate_NNt(const ElementBasis& basis) const {
        StaticMatrix<4, 4> M;
        M.setZero();
        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto p = integration_scheme_.get_point(ip);
            const auto N = shape_function(p.r, p.s);
            const auto dN_rs = shape_derivative(p.r, p.s);
            StaticVector<4> dN_da;
            StaticVector<4> dN_db;
            Precision detJ;
            Mat2 A;
            shape_gradients_physical(dN_rs, basis, dN_da, dN_db, detJ, A);
            M += p.w * detJ * (N * N.transpose());
        }
        return M;
    }

    MapMatrix mass(Precision* buffer) override {
        const ElementBasis basis = reference_basis();
        const Precision    rho   = this->get_density(false);

        Mat24 mass_matrix;
        mass_matrix.setZero();

        if (rho > Precision(0)) {
            const Precision h            = this->get_section()->thickness_;
            const Precision mass_surface = rho * h;
            const Precision inertia_rot  =
                rho * h * h * h / Precision(12);
            const auto MNN = integrate_NNt(basis);

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

    Vec3 global_point_current(Precision r, Precision s) const {
        const auto N   = shape_function(r, s);
        const auto pos = node_coords_current_6();
        Vec3 x = Vec3::Zero();
        for (Index i = 0; i < 4; ++i) {
            x += N(i) * pos.template block<1, 3>(i, 0).transpose();
        }
        return x;
    }

    Precision integrate_scalar_field(bool scale_by_density, const ScalarField& field) override {
        const ElementBasis basis = reference_basis();
        const Precision h = this->get_section()->thickness_;
        const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

        Precision result = Precision(0);
        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto p = integration_scheme_.get_point(ip);
            const auto dN_rs = shape_derivative(p.r, p.s);
            StaticVector<4> dN_da;
            StaticVector<4> dN_db;
            Precision detJ;
            Mat2 A;
            shape_gradients_physical(dN_rs, basis, dN_da, dN_db, detJ, A);
            result += field(global_point_current(p.r, p.s)) * (rho * h * p.w * detJ);
        }
        return result;
    }

    Vec3 integrate_vector_field(bool scale_by_density, const VecField& field) override {
        const ElementBasis basis = reference_basis();
        const Precision h = this->get_section()->thickness_;
        const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

        Vec3 result = Vec3::Zero();
        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto p = integration_scheme_.get_point(ip);
            const auto dN_rs = shape_derivative(p.r, p.s);
            StaticVector<4> dN_da;
            StaticVector<4> dN_db;
            Precision detJ;
            Mat2 A;
            shape_gradients_physical(dN_rs, basis, dN_da, dN_db, detJ, A);
            result += field(global_point_current(p.r, p.s)) * (rho * h * p.w * detJ);
        }
        return result;
    }

    void integrate_vector_field(Field& node_loads, bool scale_by_density, const VecField& field) override {
        const ElementBasis basis = reference_basis();
        const Precision h = this->get_section()->thickness_;
        const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto p = integration_scheme_.get_point(ip);
            const auto N = shape_function(p.r, p.s);
            const auto dN_rs = shape_derivative(p.r, p.s);
            StaticVector<4> dN_da;
            StaticVector<4> dN_db;
            Precision detJ;
            Mat2 A;
            shape_gradients_physical(dN_rs, basis, dN_da, dN_db, detJ, A);

            const Vec3 f = field(global_point_current(p.r, p.s)) * (rho * h * p.w * detJ);
            for (Index i = 0; i < 4; ++i) {
                const Index node_id = static_cast<Index>(this->node_ids[i]);
                for (Index d = 0; d < 3; ++d) {
                    node_loads(node_id, d) += N(i) * f(d);
                }
            }
        }
    }

    Mat3 integrate_tensor_field(bool scale_by_density, const TenField& field) override {
        const ElementBasis basis = reference_basis();
        const Precision h = this->get_section()->thickness_;
        const Precision rho = scale_by_density ? this->get_density(true) : Precision(1);

        Mat3 result = Mat3::Zero();
        for (Index ip = 0; ip < integration_scheme_.count(); ++ip) {
            const auto p = integration_scheme_.get_point(ip);
            const auto dN_rs = shape_derivative(p.r, p.s);
            StaticVector<4> dN_da;
            StaticVector<4> dN_db;
            Precision detJ;
            Mat2 A;
            shape_gradients_physical(dN_rs, basis, dN_da, dN_db, detJ, A);
            result += field(global_point_current(p.r, p.s)) * (rho * h * p.w * detJ);
        }
        return result;
    }

    bool compute_shell_section_forces(Field& resultants,
                                      Field& contribution_count,
                                      const Field& displacement) override {
        (void) displacement;
        const ElementBasis basis = reference_basis();
        const CurrentState state = current_state(basis);
        const Mat8 H = resultant_stiffness();

        const RowMatrix rst = stress_strain_nodal_rst();
        for (Index i = 0; i < 4; ++i) {
            Vec8 eps;
            const Vec8 res = shell_resultant_at(state, basis, H, rst(i, 0), rst(i, 1), &eps);
            const Index node_id = static_cast<Index>(this->node_ids[i]);
            for (Index a = 0; a < std::min<Index>(resultants.components, num_strains); ++a) {
                resultants(node_id, a) += res(a);
            }
            contribution_count(node_id, 0) += Precision(1);
        }
        return true;
    }
};

} // namespace fem::model
