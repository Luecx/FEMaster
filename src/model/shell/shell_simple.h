//
// Created by f_eggers on 11.12.2024.
//

#ifndef SHELL_SIMPLE_H
#define SHELL_SIMPLE_H

#include "../../material/stress/volume_stress_cauchy.h"
#include "../geometry/surface/surface.h"
#include "../geometry/surface/surface3.h"
#include "../geometry/surface/surface4.h"
#include "../geometry/surface/surface6.h"
#include "../geometry/surface/surface8.h"
#include "shell.h"

#include <algorithm>
#include <functional>

// this file is partially based on the implementation of
// https://github.com/JWock82/Pynite/blob/main/Archived/MITC4.py
//

namespace fem::model {

template<Index N, typename SFType, math::quadrature::Domain INT_D, math::quadrature::Order INT_O>
struct DefaultShellElement : public ShellElement<N> {
    SFType                 geometry;
    math::quadrature::Quadrature integration_scheme_;

    using Axes            = Mat3;
    using NodeCoords      = StaticMatrix<N, 3>;
    using LocalCoords     = StaticMatrix<N, 2>;
    using Jacobian        = StaticMatrix<2, 2>;
    using ShapeDerivative = StaticMatrix<N, 2>;
    using ShapeFunction   = StaticVector<N>;

    DefaultShellElement(ID p_elem_id, std::array<ID, N> p_node)
        : ShellElement<N>(p_elem_id, p_node)
        , geometry(p_node)
        , integration_scheme_(INT_D, INT_O) {}

    ~DefaultShellElement() override = default;

    Axes get_xyz_axes() {
        NodeCoords node_coords = this->node_coords_global();

        Vec3 n1  = node_coords.row(0);
        Vec3 n2  = node_coords.row(1);
        Vec3 n3  = node_coords.row(2);

        Vec3 n12 = n2 - n1;
        Vec3 n13 = n3 - n1;

        Vec3 x_axis = n12.normalized();
        Vec3 z_axis = n12.cross(n13).normalized();
        Vec3 y_axis = z_axis.cross(x_axis).normalized();

        StaticMatrix<3, 3> res;
        res.row(0) = x_axis.transpose();
        res.row(1) = y_axis.transpose();
        res.row(2) = z_axis.transpose();

        return res;
    }

    StaticMatrix<N * 6, N * 6> transformation(Mat3& axes) {
        StaticMatrix<N * 6, N * 6> res = StaticMatrix<N * 6, N * 6>::Zero();

        for (Index i = 0; i < 2 * N; ++i) {
            res.template block<3, 3>(3 * i, 3 * i) = axes;
        }

        return res;
    }

    LocalCoords get_xy_coords(Mat3& axes) {
        NodeCoords node_coords = this->node_coords_global();

        Vec3 x_axis = axes.row(0).transpose();
        Vec3 y_axis = axes.row(1).transpose();

        StaticMatrix<N, 2> xy_coords;

        for (Index i = 0; i < N; i++) {
            Vec3 ni         = node_coords.row(i);
            Vec3 n0         = node_coords.row(0);
            Vec3 n          = ni - n0;

            xy_coords(i, 0) = n.dot(x_axis);
            xy_coords(i, 1) = n.dot(y_axis);
        }

        return xy_coords;
    }

    Vec3 reference_point(Precision r, Precision s) {
        ShapeFunction H           = this->shape_function(r, s);
        NodeCoords    node_coords = this->node_coords_reference();

        Vec3 point = Vec3::Zero();

        for (Index i = 0; i < N; ++i) {
            point += H(i) * node_coords.row(i).transpose();
        }

        return point;
    }

    ShapeFunction shape_function(Precision r, Precision s) {
        return geometry.shape_function(r, s);
    }

    ShapeDerivative shape_derivative(Precision r, Precision s) {
        return geometry.shape_derivative(r, s);
    }

    Jacobian jacobian(StaticMatrix<N, 2>& shape_derivatives, StaticMatrix<N, 2>& xy_coords) {
        Jacobian jacobian;

        jacobian(0, 0) = shape_derivatives.col(0).dot(xy_coords.col(0));
        jacobian(1, 0) = shape_derivatives.col(0).dot(xy_coords.col(1));

        jacobian(0, 1) = shape_derivatives.col(1).dot(xy_coords.col(0));
        jacobian(1, 1) = shape_derivatives.col(1).dot(xy_coords.col(1));

        return jacobian;
    }

    Precision integrate_scalar_field(
        bool                                       scale_by_density,
        const typename ShellElement<N>::ScalarField& field
    ) override {
        Mat3 axes      = this->get_xyz_axes();
        auto xy_coords = this->get_xy_coords(axes);

        const auto&     scheme = this->integration_scheme();
        const Precision t      = this->get_section()->thickness_;

        Precision rho = Precision(1);

        if (scale_by_density) {
            rho = this->get_density(true);
        }

        auto      node_coords_glob = this->node_coords_global();
        Precision result           = Precision(0);

        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto      p = scheme.get_point(ip);
            const Precision r = p.r;
            const Precision s = p.s;
            const Precision w = p.w;

            ShapeFunction   H    = this->shape_function(r, s);
            ShapeDerivative dH   = this->shape_derivative(r, s);
            Jacobian        Jrs  = this->jacobian(dH, xy_coords);
            const Precision detJ = std::abs(Jrs.determinant());

            Vec3 x_ip = Vec3::Zero();

            for (Index i = 0; i < N; ++i) {
                x_ip += H(i) * node_coords_glob.row(i);
            }

            result += field(x_ip) * (rho * t * w * detJ);
        }

        return result;
    }

    Vec3 integrate_vector_field(
        bool                                    scale_by_density,
        const typename ShellElement<N>::VecField& field
    ) override {
        Mat3 axes      = this->get_xyz_axes();
        auto xy_coords = this->get_xy_coords(axes);

        const auto&     scheme = this->integration_scheme();
        const Precision t      = this->get_section()->thickness_;

        Precision rho = Precision(1);

        if (scale_by_density) {
            rho = this->get_density(true);
        }

        auto node_coords_glob = this->node_coords_global();
        Vec3 result           = Vec3::Zero();

        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto      p = scheme.get_point(ip);
            const Precision r = p.r;
            const Precision s = p.s;
            const Precision w = p.w;

            ShapeFunction   H    = this->shape_function(r, s);
            ShapeDerivative dH   = this->shape_derivative(r, s);
            Jacobian        Jrs  = this->jacobian(dH, xy_coords);
            const Precision detJ = std::abs(Jrs.determinant());

            Vec3 x_ip = Vec3::Zero();

            for (Index i = 0; i < N; ++i) {
                x_ip += H(i) * node_coords_glob.row(i);
            }

            result += field(x_ip) * (rho * t * w * detJ);
        }

        return result;
    }

    void integrate_vector_field(
        Field&                                  node_loads,
        bool                                    scale_by_density,
        const typename ShellElement<N>::VecField& field
    ) override {
        Mat3 axes      = this->get_xyz_axes();
        auto xy_coords = this->get_xy_coords(axes);

        const auto&     scheme = this->integration_scheme();
        const Precision t      = this->get_section()->thickness_;

        Precision rho = Precision(1);

        if (scale_by_density) {
            rho = this->get_density(true);
        }

        auto node_coords_glob = this->node_coords_global();

        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto      p = scheme.get_point(ip);
            const Precision r = p.r;
            const Precision s = p.s;
            const Precision w = p.w;

            ShapeFunction   H    = this->shape_function(r, s);
            ShapeDerivative dH   = this->shape_derivative(r, s);
            Jacobian        Jrs  = this->jacobian(dH, xy_coords);
            const Precision detJ = std::abs(Jrs.determinant());

            Vec3 x_ip = Vec3::Zero();

            for (Index i = 0; i < N; ++i) {
                x_ip += H(i) * node_coords_glob.row(i);
            }

            Vec3 f_ip = field(x_ip) * (rho * t * w * detJ);

            for (Index i = 0; i < N; ++i) {
                const ID        n_id = this->node_ids[i];
                const Precision a    = H(i);

                node_loads(n_id, 0) += a * f_ip(0);
                node_loads(n_id, 1) += a * f_ip(1);
                node_loads(n_id, 2) += a * f_ip(2);
            }
        }
    }

    Mat3 integrate_tensor_field(
        bool                                    scale_by_density,
        const typename ShellElement<N>::TenField& field
    ) override {
        Mat3 axes      = this->get_xyz_axes();
        auto xy_coords = this->get_xy_coords(axes);

        const auto&     scheme = this->integration_scheme();
        const Precision t      = this->get_section()->thickness_;

        Precision rho = Precision(1);

        if (scale_by_density) {
            rho = this->get_density(true);
        }

        auto node_coords_glob = this->node_coords_global();
        Mat3 result           = Mat3::Zero();

        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto      p = scheme.get_point(ip);
            const Precision r = p.r;
            const Precision s = p.s;
            const Precision w = p.w;

            ShapeFunction   H    = this->shape_function(r, s);
            ShapeDerivative dH   = this->shape_derivative(r, s);
            Jacobian        Jrs  = this->jacobian(dH, xy_coords);
            const Precision detJ = std::abs(Jrs.determinant());

            Vec3 x_ip = Vec3::Zero();

            for (Index i = 0; i < N; ++i) {
                x_ip += H(i) * node_coords_glob.row(i);
            }

            result += field(x_ip) * (rho * t * w * detJ);
        }

        return result;
    }

    StaticMatrix<3, N * 3> strain_disp_bending(
        ShapeDerivative& shape_der,
        Jacobian&        jacobian
    ) {
        Mat2 inv = jacobian.inverse();
        auto dH  = (shape_der * inv).transpose();

        StaticMatrix<3, N * 3> res {};
        res.setZero();

        for (Index i = 0; i < N; i++) {
            res(0, 3 * i)     = 0;
            res(0, 3 * i + 1) = 0;
            res(0, 3 * i + 2) = -dH(0, i);

            res(1, 3 * i)     = 0;
            res(1, 3 * i + 1) = dH(1, i);
            res(1, 3 * i + 2) = 0;

            res(2, 3 * i)     = 0;
            res(2, 3 * i + 1) = dH(0, i);
            res(2, 3 * i + 2) = -dH(1, i);
        }

        return res;
    }

    virtual StaticMatrix<2, 3 * N> strain_disp_shear(
        ShapeFunction&   shape_func,
        ShapeDerivative& shape_der,
        Jacobian&        jacobian
    ) {
        auto H  = shape_func;
        auto dH = (shape_der * jacobian.inverse()).transpose();

        StaticMatrix<2, 3 * N> res {};
        res.setZero();

        for (Index i = 0; i < N; i++) {
            res(0, 3 * i)     = dH(0, i);
            res(0, 3 * i + 2) = H(i);

            res(1, 3 * i)     = dH(1, i);
            res(1, 3 * i + 1) = -H(i);
        }

        return res;
    }

    virtual StaticMatrix<2, 3 * N> strain_disp_shear_at(
        Precision          r,
        Precision          s,
        const LocalCoords& xy_coords
    ) {
        ShapeFunction   H  = this->shape_function(r, s);
        ShapeDerivative dH = this->shape_derivative(r, s);
        Jacobian        J  = this->jacobian(dH, const_cast<LocalCoords&>(xy_coords));

        return this->strain_disp_shear(H, dH, J);
    }

    StaticMatrix<3, 2 * N> strain_disp_membrane(
        ShapeDerivative& shape_der,
        Jacobian&        jacobian
    ) {
        auto dH = (shape_der * jacobian.inverse()).transpose();

        StaticMatrix<3, 2 * N> res {};
        res.setZero();

        for (Index i = 0; i < N; i++) {
            res(0, 2 * i)     = dH(0, i);
            res(0, 2 * i + 1) = 0;

            res(1, 2 * i)     = 0;
            res(1, 2 * i + 1) = dH(1, i);

            res(2, 2 * i)     = dH(1, i);
            res(2, 2 * i + 1) = dH(0, i);
        }

        return res;
    }

    StaticMatrix<6 * N, 6 * N> stiffness_abd_shear(LocalCoords& xy_coords) {
        auto section = this->get_section();
        Mat3 axes        = get_xyz_axes();
        Mat3 shell_basis = axes.transpose();

        Precision topo_scale = Precision(1);

        if (this->_model_data && this->_model_data->element_stiffness_scale) {
            auto scale_field = this->_model_data->element_stiffness_scale;

            logging::error(
                scale_field->components == 1,
                "Field '", scale_field->name, "': element stiffness scale requires 1 component"
            );

            topo_scale = (*scale_field)(static_cast<Index>(this->elem_id));
        }

        std::function<StaticMatrix<5 * N, 5 * N>(Precision, Precision, Precision)> func_abd =
            [this, section, topo_scale, &xy_coords, &shell_basis](
                Precision r,
                Precision s,
                Precision /*t*/
            ) -> StaticMatrix<5 * N, 5 * N> {
                ShapeDerivative shape_der = this->shape_derivative(r, s);
                Jacobian        jac       = this->jacobian(shape_der, xy_coords);
                Precision       jac_det   = jac.determinant();

                auto Bm = this->strain_disp_membrane(shape_der, jac);
                auto Bb = this->strain_disp_bending(shape_der, jac);

                ShellGeneralizedStrain zero_strain;
                ShellStressResultants  zero_resultants;
                Mat8                   tangent;
                section->evaluate(
                    this->reference_point(r, s),
                    shell_basis,
                    zero_strain,
                    false,
                    zero_resultants,
                    tangent
                );

                const Mat6 E = topo_scale * tangent.template block<6, 6>(0, 0);

                StaticMatrix<6, 5 * N> B;
                B.setZero();

                B.template block<3, 2 * N>(0, 0)       = Bm;
                B.template block<3, 3 * N>(3, 2 * N)   = Bb;

                return B.transpose() * (E * B) * jac_det;
            };

        std::function<StaticMatrix<3 * N, 3 * N>(Precision, Precision, Precision)> func_shear =
            [this, section, topo_scale, &xy_coords, &shell_basis](
                Precision r,
                Precision s,
                Precision /*t*/
            ) -> StaticMatrix<3 * N, 3 * N> {
                auto Bs = this->strain_disp_shear_at(r, s, xy_coords);

                ShapeDerivative dH   = this->shape_derivative(r, s);
                Jacobian        J    = this->jacobian(dH, xy_coords);
                Precision       detJ = J.determinant();

                ShellGeneralizedStrain zero_strain;
                ShellStressResultants  zero_resultants;
                Mat8                   tangent;
                section->evaluate(
                    this->reference_point(r, s),
                    shell_basis,
                    zero_strain,
                    false,
                    zero_resultants,
                    tangent
                );

                const Mat2 E = topo_scale * tangent.template block<2, 2>(6, 6);

                return Bs.transpose() * (E * Bs) * detJ;
            };

        StaticMatrix<5 * N, 5 * N> stiff_abd   = this->integration_scheme().integrate(func_abd);
        StaticMatrix<3 * N, 3 * N> stiff_shear = this->integration_scheme().integrate(func_shear);

        Precision k_drill = Precision(1e32);

        for (Index i = 0; i < N; i++) {
            k_drill = std::min(k_drill, stiff_shear(3 * i + 1, 3 * i + 1));
            k_drill = std::min(k_drill, stiff_shear(3 * i + 2, 3 * i + 2));
        }

        k_drill /= Precision(1000);

        StaticMatrix<6 * N, 6 * N> res;
        res.setZero();

        int abd_dof_map[] { 0, 1, 2, 3, 4,};

        for (Index i = 0; i < 5 * N; i++) {
            for (Index j = 0; j < 5 * N; j++) {
                auto i_local_id   = i < 2 * N ? i / 2 : (i - 2 * N) / 3;
                auto j_local_id   = j < 2 * N ? j / 2 : (j - 2 * N) / 3;
                auto i_dof        = i < 2 * N ? i % 2 : 2 + (i - 2 * N) % 3;
                auto j_dof        = j < 2 * N ? j % 2 : 2 + (j - 2 * N) % 3;

                auto i_glob_index = i_local_id * 6 + abd_dof_map[i_dof];
                auto j_glob_index = j_local_id * 6 + abd_dof_map[j_dof];

                res(i_glob_index, j_glob_index) += stiff_abd(i, j);
            }
        }

        int shear_dof_map[] {2,3,4};

        for (Index i = 0; i < 3 * N; i++) {
            for (Index j = 0; j < 3 * N; j++) {
                auto i_local_id   = i / 3;
                auto j_local_id   = j / 3;
                auto i_dof        = i % 3;
                auto j_dof        = j % 3;

                auto i_glob_index = i_local_id * 6 + shear_dof_map[i_dof];
                auto j_glob_index = j_local_id * 6 + shear_dof_map[j_dof];

                res(i_glob_index, j_glob_index) += stiff_shear(i, j);
            }
        }

        for (Index i = 0; i < N; i++) {
            res(6 * i + 5, 6 * i + 5) = k_drill;
        }

        return res;
    }

    const math::quadrature::Quadrature& integration_scheme() const override {
        return integration_scheme_;
    }

    MapMatrix stiffness(Precision* buffer) override {
        auto axes      = get_xyz_axes();
        auto xy_coords = get_xy_coords(axes);

        MapMatrix mapped {buffer, 6 * N, 6 * N};

        auto trans = transformation(axes);
        auto stiff = stiffness_abd_shear(xy_coords);

        mapped = trans.transpose() * stiff * trans;

        return mapped;
    }

    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override {
        auto axes      = get_xyz_axes();
        auto xy_coords = get_xy_coords(axes);

        Index ip_counter = 0;

        std::function<StaticMatrix<3 * N, 3 * N>(Precision, Precision, Precision)> func_geo =
            [this, &ip_stress, ip_start_idx, &ip_counter, &xy_coords](
                Precision r,
                Precision s,
                Precision /*t*/
            ) -> StaticMatrix<3 * N, 3 * N> {
                ShapeDerivative dH_rs = this->shape_derivative(r, s);
                Jacobian        J     = this->jacobian(dH_rs, const_cast<LocalCoords&>(xy_coords));
                Mat2            Jinv  = J.inverse();
                Precision       detJ  = J.determinant();

                auto dH_xy = (dH_rs * Jinv).transpose();

                const Index ip_row = static_cast<Index>(ip_start_idx) + ip_counter++;
                const Vec6  v      = ip_stress.row_vec6(ip_row);

                const Precision nxx = v(0);
                const Precision nyy = v(1);
                const Precision nxy = v(5);

                Eigen::Matrix<Precision, 2, 2> Nmat;
                Nmat << nxx, nxy,
                        nxy, nyy;

                StaticMatrix<3 * N, 3 * N> Kloc;
                Kloc.setZero();

                for (Index a = 0; a < N; ++a) {
                    Eigen::Matrix<Precision, 2, 1> ga;
                    ga << dH_xy(0, a), dH_xy(1, a);

                    for (Index b = 0; b < N; ++b) {
                        Eigen::Matrix<Precision, 2, 1> gb;
                        gb << dH_xy(0, b), dH_xy(1, b);

                        const Precision kab = (ga.transpose() * Nmat * gb)(0, 0) * detJ;

                        for (Index c = 0; c < 3; ++c) {
                            Kloc(3 * a + c, 3 * b + c) += kab;
                        }
                    }
                }

                return Kloc;
            };

        StaticMatrix<3 * N, 3 * N> Kg_trans_local = this->integration_scheme().integrate(func_geo);
        Kg_trans_local = Precision(0.5) * (Kg_trans_local + Kg_trans_local.transpose());

        StaticMatrix<6 * N, 6 * N> Kg6;
        Kg6.setZero();

        int dof_map[] {0, 1, 2};

        for (Index i = 0; i < 3 * N; ++i) {
            for (Index j = 0; j < 3 * N; ++j) {
                const Index i_node = i / 3;
                const Index i_ldof = i % 3;
                const Index j_node = j / 3;
                const Index j_ldof = j % 3;

                const Index I = 6 * i_node + dof_map[i_ldof];
                const Index J = 6 * j_node + dof_map[j_ldof];

                Kg6(I, J) = Kg_trans_local(i, j);
            }
        }

        MapMatrix mapped {buffer, 6 * N, 6 * N};

        auto T = transformation(axes);
        mapped = T.transpose() * Kg6 * T;

        return mapped;
    }

    StaticMatrix<N, N> integrate_NNt(const LocalCoords& xy_coords) {
        std::function<StaticMatrix<N, N>(Precision, Precision, Precision)> fn =
            [this, &xy_coords](
                Precision r,
                Precision s,
                Precision /*t*/
            ) -> StaticMatrix<N, N> {
                ShapeFunction   H    = this->shape_function(r, s);
                ShapeDerivative dH   = this->shape_derivative(r, s);
                Jacobian        J    = this->jacobian(dH, const_cast<LocalCoords&>(xy_coords));
                Precision       detJ = J.determinant();

                return (H * H.transpose()) * detJ;
            };

        return this->integration_scheme().integrate(fn);
    }

    MapMatrix mass(Precision* buffer) override {
        auto axes      = get_xyz_axes();
        auto xy_coords = get_xy_coords(axes);

        const Precision rho = this->get_density(false);

        if (rho <= Precision(0)) {
            MapMatrix mapped(buffer, 6 * N, 6 * N);
            mapped.setZero();
            return mapped;
        }

        const Precision h  = this->get_section()->thickness_;
        const Precision mt = rho * h;
        const Precision mr = rho * (h * h * h / Precision(12));

        StaticMatrix<N, N> M_NN = integrate_NNt(xy_coords);

        StaticMatrix<6 * N, 6 * N> Mloc;
        Mloc.setZero();

        for (Index i = 0; i < N; ++i) {
            for (Index j = 0; j < N; ++j) {
                const Precision mij = mt * M_NN(i, j);

                for (int k = 0; k < 3; ++k) {
                    const Index I = 6 * i + k;
                    const Index J = 6 * j + k;

                    Mloc(I, J) += mij;
                }
            }
        }

        for (Index i = 0; i < N; ++i) {
            for (Index j = 0; j < N; ++j) {
                const Precision mij = mr * M_NN(i, j);

                for (int k = 0; k < 2; ++k) {
                    const Index I = 6 * i + 3 + k;
                    const Index J = 6 * j + 3 + k;

                    Mloc(I, J) += mij;
                }
            }
        }

        {
            Precision avg_m_node = Precision(0);

            for (Index i = 0; i < N; ++i) {
                avg_m_node += mt * M_NN(i, i);
            }

            avg_m_node /= std::max<Precision>(Precision(1), static_cast<Precision>(N));

            const Precision m_drill = Precision(1e-6) * avg_m_node;

            for (Index i = 0; i < N; ++i) {
                const Index I = 6 * i + 5;
                Mloc(I, I) += m_drill;
            }
        }

        MapMatrix mapped(buffer, 6 * N, 6 * N);

        auto T = transformation(axes);
        mapped = T.transpose() * Mloc * T;

        return mapped;
    }

    Precision volume() override {
        const Precision h = this->get_section()->thickness_;

        logging::error(this->_model_data != nullptr, "no model data assigned to element ", this->elem_id);
        logging::error(this->_model_data->positions != nullptr, "positions field not set in model data");

        const auto& node_coords_system = *this->_model_data->positions;
        const Precision A              = geometry.area(node_coords_system);

        return h * A;
    }

    virtual VolumeStressCauchy compute_stress_at(Field& displacement, Vec3& xyz) {
        // local coordinates extracted
        const Precision r = xyz(0);
        const Precision s = xyz(1);
        const Precision t = xyz(2);

        // axes of the current element
        Mat3            axes        = get_xyz_axes();
        Mat3            shell_basis = axes.transpose();
        LocalCoords     xy_coords   = get_xy_coords(axes);
        ShapeDerivative shape_der   = this->shape_derivative(r, s);
        Jacobian        jac         = this->jacobian(shape_der, xy_coords);

        // extract thickness
        const Precision h = this->get_section()->thickness_;

        // for ever node, we need to store the in-plane dispalacement (u,v), bending components (w,phi1, phi2) and the
        // shear components (w,phi1,phi2)
        StaticVector<2 * N> disp_membrane;
        StaticVector<3 * N> disp_bending;
        StaticVector<3 * N> disp_shear;

        disp_membrane.setZero();
        disp_bending .setZero();
        disp_shear   .setZero();

        for (Index i = 0; i < N; i++) {
            ID node_id = this->nodes()[i];

            Vec6 displacement_glob = displacement.row_vec6(static_cast<Index>(node_id));
            Vec3 disp_xyz          = displacement_glob.head(3);
            Vec3 disp_rot          = displacement_glob.tail(3);

            disp_xyz = axes * disp_xyz;
            disp_rot = axes * disp_rot;

            // membrane part
            disp_membrane(2 * i)     = disp_xyz(0);
            disp_membrane(2 * i + 1) = disp_xyz(1);

            // bending dofs
            disp_bending(3 * i)      = disp_xyz(2);
            disp_bending(3 * i + 1)  = disp_rot(0);
            disp_bending(3 * i + 2)  = disp_rot(1);

            // shear dofs
            disp_shear(3 * i)        = disp_xyz(2);
            disp_shear(3 * i + 1)    = disp_rot(0);
            disp_shear(3 * i + 2)    = disp_rot(1);
        }

        // in case of topology optimization, this elements stiffness may be scaled
        Precision topo_scale = Precision(1);
        if (this->_model_data && this->_model_data->element_stiffness_scale) {
            auto scale_field = this->_model_data->element_stiffness_scale;
            topo_scale = (*scale_field)(static_cast<Index>(this->elem_id));
        }

        // get the B matrices for membrane, bending and shear
        auto B_membrane = this->strain_disp_membrane(shape_der, jac);
        auto B_bending  = this->strain_disp_bending(shape_der, jac);
        auto B_shear    = this->strain_disp_shear_at(r, s, xy_coords);

        // compute the actual generalised strain values
        Vec3 eps_element   = B_membrane * disp_membrane;
        Vec3 kappa_element = B_bending  * disp_bending;
        Vec2 gamma_element = B_shear    * disp_shear;

        // store them inside a ShellGeneralizedStrain object
        constexpr Index membrane_start =
            static_cast<Index>(ShellGeneralizedStrain::Component::EpsilonXX);
        constexpr Index curvature_start =
            static_cast<Index>(ShellGeneralizedStrain::Component::KappaXX);
        constexpr Index shear_start =
            static_cast<Index>(ShellGeneralizedStrain::Component::GammaXZ);

        ShellGeneralizedStrain generalized_strain;
        generalized_strain.values().template segment<3>(membrane_start)  = eps_element;
        generalized_strain.values().template segment<3>(curvature_start) = kappa_element;
        generalized_strain.values().template segment<2>(shear_start)     = gamma_element;

        const Precision z = t * h / Precision(2);

        VolumeStressCauchy stress = this->get_section()->evaluate_output_stress(
            reference_point(r, s),
            shell_basis,
            generalized_strain,
            z,
            false
        );
        stress.voigt() *= topo_scale;

        return stress;
    }

    bool compute_shell_section_forces(Field& resultants,
                                      Field& contribution_count,
                                      const Field& displacement) override {
        // get local axes
        Mat3 axes        = get_xyz_axes();
        Mat3 shell_basis = axes.transpose();

        Precision topo_scale = Precision(1);

        if (this->_model_data && this->_model_data->element_stiffness_scale) {
            auto scale_field = this->_model_data->element_stiffness_scale;
            topo_scale       = (*scale_field)(static_cast<Index>(this->elem_id));
        }

        StaticVector<2 * N> disp_membrane;
        StaticVector<3 * N> disp_bending;
        StaticVector<3 * N> disp_shear;

        disp_membrane.setZero();
        disp_bending.setZero();
        disp_shear.setZero();

        for (Index a = 0; a < N; ++a) {
            ID node_id = this->nodes()[a];

            Vec6 displacement_glob = displacement.row_vec6(static_cast<Index>(node_id));
            Vec3 disp_xyz          = displacement_glob.head(3);
            Vec3 disp_rot          = displacement_glob.tail(3);

            disp_xyz = axes * disp_xyz;
            disp_rot = axes * disp_rot;

            disp_membrane(2 * a)     = disp_xyz(0);
            disp_membrane(2 * a + 1) = disp_xyz(1);

            disp_bending(3 * a)      = disp_xyz(2);
            disp_bending(3 * a + 1)  = disp_rot(0);
            disp_bending(3 * a + 2)  = disp_rot(1);

            disp_shear(3 * a)        = disp_xyz(2);
            disp_shear(3 * a + 1)    = disp_rot(0);
            disp_shear(3 * a + 2)    = disp_rot(1);
        }

        StaticMatrix<N, 2> coords    = this->geometry.node_coords_local();
        LocalCoords        xy_coords = get_xy_coords(axes);

        for (Index i = 0; i < N; ++i) {
            ID node_id = this->nodes()[i];

            Precision r = coords(i, 0);
            Precision s = coords(i, 1);

            ShapeDerivative shape_der = this->shape_derivative(r, s);
            Jacobian        jac       = this->jacobian(shape_der, xy_coords);

            auto B_membrane = this->strain_disp_membrane(shape_der, jac);
            auto B_bending  = this->strain_disp_bending(shape_der, jac);
            auto B_shear    = this->strain_disp_shear_at(r, s, xy_coords);

            Vec3 eps_element   = B_membrane * disp_membrane;
            Vec3 kappa_element = B_bending  * disp_bending;
            Vec2 gamma_element = B_shear    * disp_shear;

            constexpr Index membrane_start =
                static_cast<Index>(ShellGeneralizedStrain::Component::EpsilonXX);
            constexpr Index curvature_start =
                static_cast<Index>(ShellGeneralizedStrain::Component::KappaXX);
            constexpr Index shear_start =
                static_cast<Index>(ShellGeneralizedStrain::Component::GammaXZ);

            ShellGeneralizedStrain generalized_strain;
            generalized_strain.values().template segment<3>(membrane_start) = eps_element;
            generalized_strain.values().template segment<3>(curvature_start) = kappa_element;
            generalized_strain.values().template segment<2>(shear_start) = gamma_element;
            const ShellStressResultants output_resultants = this->get_section()->evaluate_output_resultants(
                reference_point(r, s),
                shell_basis,
                generalized_strain,
                false
            );

            const Vec8 values = topo_scale * output_resultants.values();

            const Index node_idx = static_cast<Index>(node_id);

            for (Index component = 0; component < 8; ++component) {
                resultants(node_idx, component) += values(component);
            }
            contribution_count(node_idx, 0) += Precision(1);
        }
        return true;
    }

    RowMatrix stress_strain_nodal_rst() override {
        StaticMatrix<N, 2> coords = this->geometry.node_coords_local();
        RowMatrix rst(N, 3);
        rst.setZero();
        for (Index i = 0; i < N; ++i) {
            rst(i, 0) = coords(i, 0);
            rst(i, 1) = coords(i, 1);
            rst(i, 2) = Precision(0);
        }
        return rst;
    }

    RowMatrix stress_strain_ip_rst() override {
        const auto& scheme = this->integration_scheme();
        RowMatrix rst(scheme.count(), 3);
        rst.setZero();
        for (Index i = 0; i < scheme.count(); ++i) {
            rst(i, 0) = scheme.get_point(i).r;
            rst(i, 1) = scheme.get_point(i).s;
            rst(i, 2) = Precision(0);
        }
        return rst;
    }

    void compute_stress_strain(Field* strain,
                               Field* stress_out,
                               const Field& displacement,
                               const RowMatrix& rst,
                               int offset,
                               bool use_green_lagrange_nl) override {
        logging::error(!use_green_lagrange_nl,
            "ShellElement: nonlinear stress/strain evaluation is not implemented yet for element ", this->elem_id);
        logging::error(strain != nullptr || stress_out != nullptr,
            "ShellElement: compute_stress_strain requires at least one output field");
        logging::error(rst.cols() >= 3,
            "ShellElement: stress/strain evaluation coordinates require at least 3 columns");

        Field& displacement_mut = const_cast<Field&>(displacement);
        for (Eigen::Index i = 0; i < rst.rows(); ++i) {
            Vec3 local;
            local << rst(i, 0), rst(i, 1), rst(i, 2);
            const VolumeStressCauchy stress_value = this->compute_stress_at(displacement_mut, local);
            const Index row = static_cast<Index>(offset + i);

            if (strain) {
                for (Index j = 0; j < strain->components; ++j) {
                    (*strain)(row, j) = Precision(0);
                }
            }
            if (stress_out) {
                for (Index j = 0; j < stress_out->components; ++j) {
                    (*stress_out)(row, j) = stress_value.voigt()(j);
                }
            }
        }
    }

    void compute_stress_state(Field& stress_state,
                              const Field& displacement,
                              int offset,
                              bool use_green_lagrange_nl) override {
        logging::error(!use_green_lagrange_nl,
            "ShellElement: nonlinear stress state is not implemented yet for element ", this->elem_id);
        compute_membrane_stress_state(stress_state, displacement, offset);
    }

    void compute_membrane_stress_state(
        Field& ip_stress,
        const Field& displacement,
        int    ip_offset
    ) {
        Mat3 axes        = get_xyz_axes();
        Mat3 shell_basis = axes.transpose();
        auto xy_coords   = get_xy_coords(axes);

        Precision topo_scale = Precision(1);

        if (this->_model_data && this->_model_data->element_stiffness_scale) {
            auto scale_field = this->_model_data->element_stiffness_scale;

            logging::error(
                scale_field->components == 1,
                "Field '", scale_field->name, "': element stiffness scale requires 1 component"
            );

            topo_scale = (*scale_field)(static_cast<Index>(this->elem_id));
        }

        StaticVector<2 * N> u_mem;
        u_mem.setZero();

        for (Index i = 0; i < N; ++i) {
            ID   node_id = this->nodes()[i];
            Vec6 u_glob  = displacement.row_vec6(static_cast<Index>(node_id));
            Vec3 t_glob  = u_glob.head<3>();
            Vec3 t_loc   = axes * t_glob;

            u_mem(2 * i)     = t_loc(0);
            u_mem(2 * i + 1) = t_loc(1);
        }

        const auto& scheme = this->integration_scheme();

        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const Precision r = scheme.get_point(ip).r;
            const Precision s = scheme.get_point(ip).s;

            ShapeDerivative dH_rs = this->shape_derivative(r, s);
            Jacobian        J     = this->jacobian(dH_rs, xy_coords);

            StaticMatrix<3, 2 * N> Bm = this->strain_disp_membrane(dH_rs, J);

            Vec3 eps_element = Bm * u_mem;

            constexpr Index membrane_start =
                static_cast<Index>(ShellGeneralizedStrain::Component::EpsilonXX);

            ShellGeneralizedStrain generalized_strain;
            generalized_strain.values().template segment<3>(membrane_start) = eps_element;
            ShellStressResultants  generalized_resultants;
            Mat8                   tangent;
            this->get_section()->evaluate(
                reference_point(r, s),shell_basis,generalized_strain,false,generalized_resultants,tangent
            );

            const Vec3 membrane_force_element = topo_scale * generalized_resultants.membrane();

            const Precision nxx = membrane_force_element(0);
            const Precision nyy = membrane_force_element(1);
            const Precision nxy = membrane_force_element(2);

            const Index row = static_cast<Index>(ip_offset) + ip;

            ip_stress(row, 0) = nxx;
            ip_stress(row, 1) = nyy;
            ip_stress(row, 2) = Precision(0);
            ip_stress(row, 3) = Precision(0);
            ip_stress(row, 4) = Precision(0);
            ip_stress(row, 5) = nxy;
        }
    }
};

} // namespace fem::model

#endif // SHELL_SIMPLE_H
