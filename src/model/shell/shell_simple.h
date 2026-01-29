//
// Created by f_eggers on 11.12.2024.
//

#ifndef SHELL_SIMPLE_H
#define SHELL_SIMPLE_H

#include "../../material/stress.h"
#include "../../material/isotropic_elasticity.h"
#include "../geometry/surface/surface.h"
#include "../geometry/surface/surface3.h"
#include "../geometry/surface/surface4.h"
#include "../geometry/surface/surface6.h"
#include "../geometry/surface/surface8.h"
#include "shell.h"

// this file is partially based on the implementation of
// https://github.com/JWock82/Pynite/blob/main/Archived/MITC4.py
//

namespace fem::model {

template<Index N, typename SFType, quadrature::Domain INT_D, quadrature::Order INT_O>
struct DefaultShellElement : public ShellElement<N> {
    SFType                 geometry;
    quadrature::Quadrature integration_scheme_;

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

    // when assuming a planar element, the normal is constant
    // we can also project the shell and use a local x,y system
    // this must not be confused with the local r,s system which is used for integration
    Axes get_xyz_axes() {
        NodeCoords node_coords = this->node_coords_global();

        Vec3       n1          = node_coords.row(0);
        Vec3       n2          = node_coords.row(1);
        Vec3       n3          = node_coords.row(2);

        Vec3       n12         = n2 - n1;
        Vec3       n13         = n3 - n1;

        Vec3       x_axis      = n12 / n12.norm();
        Vec3       z_axis      = n12.cross(n13) / (n12.cross(n13)).norm();
        Vec3       y_axis      = z_axis.cross(x_axis);

        // for each node calculate the local x,y coordinates by projecting onto the axes
        StaticMatrix<3, 3> res;
        res.row(0) = x_axis.transpose();
        res.row(1) = y_axis.transpose();
        res.row(2) = z_axis.transpose();
        return res;
    }

    // Assumption is that the element is planar
    // deviations of that will produce inaccurate results
    // finding the transformation matrix is straight forward
    // the x-axis will be aligned with the direction from node 1 to node 2
    // the z-axis will be perpendicular to the plane of the element
    // the y-axis will be the cross product of the other two
    StaticMatrix<N * 6, N * 6> transformation(Mat3& axes) {
        StaticMatrix<N * 6, N * 6> res = StaticMatrix<N * 6, N * 6>::Zero();
        for (int i = 0; i < 2 * N; ++i) {    // Loop over the 8 nodes
            res.template block<3, 3>(3 * i, 3 * i) = axes;
        }
        return res;
    }

    LocalCoords get_xy_coords(Mat3& axes) {
        NodeCoords node_coords = this->node_coords_global();

        Vec3       x_axis      = axes.row(0).transpose();
        Vec3       y_axis      = axes.row(1).transpose();

        // for each node calculate the local x,y coordinates by projecting onto the axes
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

    ShapeFunction shape_function(Precision r, Precision s) {
        return geometry.shape_function(r, s);
    }

    ShapeDerivative shape_derivative(Precision r, Precision s) {
        return geometry.shape_derivative(r, s);
    }

    // function which relates derivatives in the local r,s to the local x,y system
    // the entries will contain the mapping
    // [dx/dr, dx/ds]
    // [dy/dr, dy/ds]
    // x = N1*x1 + N2*x2 + N3*x3 + N4*x4
    // y = N1*y1 + N2*y2 + N3*y3 + N4*y4
    // hence the derivatives are
    // dx/dr = dN1/dr*x1 + dN2/dr*x2 + dN3/dr*x3 + dN4/dr*x4
    Jacobian jacobian(StaticMatrix<N, 2>& shape_derivatives, StaticMatrix<N, 2>& xy_coords) {
        // computing the jacobian is done by computing the derivatives of the shape functions
        // with the local x,y system
        Jacobian jacobian;
        jacobian(0, 0) = shape_derivatives.col(0).dot(xy_coords.col(0));    // dx/dr
        jacobian(1, 0) = shape_derivatives.col(0).dot(xy_coords.col(1));    // dy/dr

        jacobian(0, 1) = shape_derivatives.col(1).dot(xy_coords.col(0));    // dx/ds
        jacobian(1, 1) = shape_derivatives.col(1).dot(xy_coords.col(1));    // dy/ds
        return jacobian;
    }

    StaticMatrix<3, N * 3> strain_disp_bending(ShapeDerivative& shape_der, Jacobian& jacobian) {
        Mat2                   inv = jacobian.inverse();
        auto                   dH  = (shape_der * inv).transpose();

        StaticMatrix<3, N * 3> res {};

        // dofs are displacement in z, rotation around x, rotation around y
        // its only a function of the rotational dofs (second derivative of displacement)
        for (int i = 0; i < N; i++) {
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
        // res <<
        //     0,    0,     -dH(0, 0), 0,    0,     -dH(0, 1), 0,    0,     -dH(0, 2), 0,    0,     -dH(0, 3),
        //     0, dH(1, 0),     0,     0, dH(1, 1),     0,     0, dH(1, 2),     0,     0, dH(1, 3),     0    ,
        //     0, dH(0, 0), -dH(1, 0), 0, dH(0, 1), -dH(1, 1), 0, dH(0, 2), -dH(1, 2), 0, dH(0, 3), -dH(1, 3);
        return res;
    }

    virtual StaticMatrix<2, 3 * N>
        strain_disp_shear(ShapeFunction& shape_func, ShapeDerivative& shape_der, Jacobian& jacobian) {
        auto H  = shape_func;
        auto dH = (shape_der * jacobian.inverse()).transpose();

        // StaticMatrix<2, 12> res{};
        // // dofs are displacement in z, rotation around x, rotation around y
        // // its only a function of the rotational dofs (second derivative of displacement)
        // res <<
        //     dH(0,0),  0, H(0), dH(0,1),  0, H(1), dH(0,2),  0, H(2), dH(0,3),  0, H(3),
        //     dH(1,0),-H(0),  0, dH(1,1),-H(1),  0, dH(1,2),-H(2),  0, dH(1,3),-H(3), 0;
        // return res;
        StaticMatrix<2, 3 * N> res {};
        for (int i = 0; i < N; i++) {
            res(0, 3 * i)     = dH(0, i);
            res(0, 3 * i + 1) = 0;
            res(0, 3 * i + 2) = H(i);

            res(1, 3 * i)     = dH(1, i);
            res(1, 3 * i + 1) = -H(i);
            res(1, 3 * i + 2) = 0;
        }
        return res;
    }

    virtual StaticMatrix<2, 3 * N> strain_disp_shear_at(Precision r, Precision s, const LocalCoords& xy_coords) {
        // Default = bisherige (nicht-MITC) Variante
        ShapeFunction   H  = this->shape_function(r, s);
        ShapeDerivative dH = this->shape_derivative(r, s);
        Jacobian        J  = this->jacobian(dH, const_cast<LocalCoords&>(xy_coords));
        return this->strain_disp_shear(H, dH, J);
    }

    StaticMatrix<3, 2 * N> strain_disp_membrane(ShapeDerivative& shape_der, Jacobian& jacobian) {
        auto dH = (shape_der * jacobian.inverse()).transpose();
        // StaticMatrix<3, 8> res{};
        // res << dH(0,0), 0      , dH(0,1), 0      , dH(0,2), 0      , dH(0,3), 0      ,
        //            0  , dH(1,0), 0      , dH(1,1), 0      , dH(1,2), 0      , dH(1,3),
        //        dH(1,0), dH(0,0), dH(1,1), dH(0,1), dH(1,2), dH(0,2), dH(1,3), dH(0,3);
        StaticMatrix<3, 2 * N> res {};
        for (int i = 0; i < N; i++) {
            res(0, 2 * i)     = dH(0, i);
            res(0, 2 * i + 1) = 0;

            res(1, 2 * i)     = 0;
            res(1, 2 * i + 1) = dH(1, i);

            res(2, 2 * i)     = dH(1, i);
            res(2, 2 * i + 1) = dH(0, i);
        }
        return res;
    }

    StaticMatrix<6 * N, 6 * N> stiffness_bending(LocalCoords& xy_coords) {

        auto mat_bend  = this->get_elasticity()->get_bend(this->get_section()->thickness);
        auto mat_shear = this->get_elasticity()->get_shear(this->get_section()->thickness);

        Precision topo_scale = Precision(1);
        if (this->_model_data && this->_model_data->element_stiffness_scale) {
            auto scale_field = this->_model_data->element_stiffness_scale;
            logging::error(scale_field->components == 1,
                           "Field '", scale_field->name, "': element stiffness scale requires 1 component");
            topo_scale = (*scale_field)(static_cast<Index>(this->elem_id));
        }
        mat_bend *= topo_scale;
        mat_shear *= topo_scale;

        std::function<StaticMatrix<3 * N, 3 * N>(Precision, Precision, Precision)> func_bend =
            [this, &mat_bend, &xy_coords](Precision r, Precision s, Precision t) -> StaticMatrix<3 * N, 3 * N> {
            ShapeDerivative shape_der = this->shape_derivative(r, s);
            Jacobian        jac       = this->jacobian(shape_der, xy_coords);
            Precision       jac_det   = jac.determinant();
            auto            B         = this->strain_disp_bending(shape_der, jac);
            auto            E         = mat_bend;
            return B.transpose() * (E * B) * jac_det;
        };

        std::function<StaticMatrix<3 * N, 3 * N>(Precision, Precision, Precision)> func_shear =
            [this, &mat_shear, &xy_coords](Precision r, Precision s, Precision t) -> StaticMatrix<3 * N, 3 * N> {
            Precision jac_det;                                             // nur fürs Interface hier nicht mehr nötig
            auto      Bs = this->strain_disp_shear_at(r, s, xy_coords);    // <-- NEU
            // wir brauchen dennoch det(J) fürs dA:
            ShapeDerivative dH   = this->shape_derivative(r, s);
            Jacobian        J    = this->jacobian(dH, xy_coords);
            Precision       detJ = J.determinant();
            auto            E    = mat_shear;
            return Bs.transpose() * (E * Bs) * detJ;
        };

        StaticMatrix<3 * N, 3 * N> stiff_bend  = this->integration_scheme().integrate(func_bend);
        StaticMatrix<3 * N, 3 * N> stiff_shear = this->integration_scheme().integrate(func_shear);
        StaticMatrix<3 * N, 3 * N> stiff_local = stiff_bend + stiff_shear;

        Precision                  k_drill     = 1e32;
        for (int i = 0; i < N; i++) {
            k_drill = std::min(k_drill, stiff_local(3 * i + 1, 3 * i + 1));
            k_drill = std::min(k_drill, stiff_local(3 * i + 2, 3 * i + 2));
        }
        k_drill /= 1000;

        StaticMatrix<6 * N, 6 * N> res;
        res.setZero();

        int dof_map[] {
            2,    // z
            3,    // rx
            4,    // ry
        };

        for (Index i = 0; i < 3 * N; i++) {
            for (Index j = 0; j < 3 * N; j++) {
                auto i_local_id                 = i / 3;
                auto j_local_id                 = j / 3;
                auto i_dof                      = i % 3;
                auto j_dof                      = j % 3;

                auto i_glob_index               = i_local_id * 6 + dof_map[i_dof];
                auto j_glob_index               = j_local_id * 6 + dof_map[j_dof];

                res(i_glob_index, j_glob_index) = stiff_local(i, j);
            }
        }

        for (int i = 0; i < N; i++) {
            res(6 * i + 5, 6 * i + 5) = k_drill;
        }

        return res;
    }

    StaticMatrix<6 * N, 6 * N> stiffness_membrane(LocalCoords& xy_coords) {
        auto mat_membrane = this->get_elasticity()->get_memb();

        Precision topo_scale = Precision(1);
        if (this->_model_data && this->_model_data->element_stiffness_scale) {
            auto scale_field = this->_model_data->element_stiffness_scale;
            logging::error(scale_field->components == 1,
                           "Field '", scale_field->name, "': element stiffness scale requires 1 component");
            topo_scale = (*scale_field)(static_cast<Index>(this->elem_id));
        }
		mat_membrane *= topo_scale;

        std::function<StaticMatrix<2 * N, 2 * N>(Precision, Precision, Precision)> func_membrane =
            [this, &mat_membrane, &xy_coords](Precision r, Precision s, Precision t) -> StaticMatrix<2 * N, 2 * N> {
            ShapeDerivative shape_der = this->shape_derivative(r, s);
            Jacobian        jac       = this->jacobian(shape_der, xy_coords);
            Precision       jac_det   = jac.determinant();
            auto            B         = this->strain_disp_membrane(shape_der, jac);
            auto            E         = mat_membrane;
            return B.transpose() * (E * B) * jac_det;
        };

        StaticMatrix<2 * N, 2 * N> stiff_membrane = this->integration_scheme().integrate(func_membrane);
        stiff_membrane *= this->get_section()->thickness;

        StaticMatrix<6 * N, 6 * N> res;
        res.setZero();

        int dof_map[] {0, 1};
        for (Index i = 0; i < 2 * N; i++) {
            for (Index j = 0; j < 2 * N; j++) {
                auto i_local_id                 = i / 2;
                auto j_local_id                 = j / 2;
                auto i_dof                      = i % 2;
                auto j_dof                      = j % 2;

                auto i_glob_index               = i_local_id * 6 + dof_map[i_dof];
                auto j_glob_index               = j_local_id * 6 + dof_map[j_dof];

                res(i_glob_index, j_glob_index) = stiff_membrane(i, j);
            }
        }
        return res;
    }

    const quadrature::Quadrature& integration_scheme() const override {
        return integration_scheme_;
    }

    MapMatrix stiffness(Precision* buffer) override {
        // compute axes and local coordinates
        auto      axes      = get_xyz_axes();
        auto      xy_coords = get_xy_coords(axes);

        MapMatrix mapped {buffer, 6 * N, 6 * N};
        auto      trans = transformation(axes);

        auto stiff = stiffness_bending(xy_coords) + stiffness_membrane(xy_coords);
        mapped     = trans.transpose() * stiff * trans;
        return mapped;
    }

    MapMatrix stiffness_geom(Precision* buffer, const Field& ip_stress, int ip_start_idx) override {
        // 1) lokale Lamina-Achsen & xy-Koordinaten (wie bei Steifigkeit/Masse)
        auto axes      = get_xyz_axes();
        auto xy_coords = get_xy_coords(axes);

        // 2) Integrand: ∫ G^T Nmat G detJ dA, nur auf (w)-DOFs im [w,rx,ry]-Block
        Index                                                                      ip_counter = 0;

        std::function<StaticMatrix<3 * N, 3 * N>(Precision, Precision, Precision)> func_geo =
            [this, &ip_stress, ip_start_idx, &ip_counter, &xy_coords](Precision r,
                                                                      Precision s,
                                                                      Precision /*t*/) -> StaticMatrix<3 * N, 3 * N> {
            // --- Jacobian & dN/dx,y ---
            ShapeDerivative dH_rs = this->shape_derivative(r, s);    // (N×2): dN/dr, dN/ds
            Jacobian        J     = this->jacobian(dH_rs, const_cast<LocalCoords&>(xy_coords));
            Mat2            Jinv  = J.inverse();
            Precision       detJ  = J.determinant();
            auto            dH_xy = (dH_rs * Jinv).transpose();    // (2×N): [dN/dx; dN/dy]

            // --- Membran-Resultierende am IP (aus ip_stress) ---
            // ip_stress: [xx, yy, zz, yz, zx, xy] (wie bei Solids)
            const Index ip_row = static_cast<Index>(ip_start_idx) + ip_counter++;
            const Vec6 v = ip_stress.row_vec6(ip_row);

            const Precision h = this->get_section()->thickness;
            // Falls ip_stress bereits Resultierende enthält, setze scale=1
            const Precision scale = 1.0;    // oder 1.0, je nachdem was du speicherst
            const Precision nxx   = scale * v(0);
            const Precision nyy   = scale * v(1);
            const Precision nxy   = scale * v(5);

            // Nmat (2x2)
            Eigen::Matrix<Precision, 2, 2> Nmat;
            Nmat << nxx, nxy, nxy, nyy;

            // --- G-Operator (2 × 3N) auf [w,rx,ry]-Blöcke:
            // Zeile 0 mappt ∂w/∂x, Zeile 1 mappt ∂w/∂y; nur die w-Spalten (Index 3*i+0) sind belegt.
            StaticMatrix<2, 3 * N> G;
            G.setZero();
            for (Index i = 0; i < N; ++i) {
                const Index c_w = 3 * i + 0;
                G(0, c_w)       = dH_xy(0, i);    // dw/dx
                G(1, c_w)       = dH_xy(1, i);    // dw/dy
            }

            // --- Lokaler 3N×3N-Block (nur w–w koppelt) ---
            StaticMatrix<3 * N, 3 * N> Kloc = G.transpose() * (Nmat * G) * detJ;
            return Kloc;
        };

        // 3) Integrieren (nur der 3N×3N [w,rx,ry]-Teil)
        StaticMatrix<3 * N, 3 * N> Kg_bend_local = this->integration_scheme().integrate(func_geo);

        // 4) In 6N-Abbildung (z, rx, ry) einsortieren – identisch zu deiner Biege/Schub-Assembly
        StaticMatrix<6 * N, 6 * N> Kg6;
        Kg6.setZero();
        int dof_map[] {2, 3, 4};    // z, rx, ry
        for (Index i = 0; i < 3 * N; ++i) {
            for (Index j = 0; j < 3 * N; ++j) {
                const Index i_node = i / 3, i_ldof = i % 3;
                const Index j_node = j / 3, j_ldof = j % 3;
                const Index I = 6 * i_node + dof_map[i_ldof];
                const Index J = 6 * j_node + dof_map[j_ldof];
                Kg6(I, J)     = Kg_bend_local(i, j);
            }
        }

        // 5) (optional) kein Drill-Anteil in Kg

        // 6) In globale Koordinaten drehen (wie bei 'stiffness')
        MapMatrix mapped {buffer, 6 * N, 6 * N};
        auto      T = transformation(axes);    // 6N×6N
        mapped      = T.transpose() * Kg6 * T;
        return mapped;
    }

    // Hilfsfunktion: ∫_A N N^T det(J) dA  (N×N)
    StaticMatrix<N, N> integrate_NNt(const LocalCoords& xy_coords) {
        std::function<StaticMatrix<N, N>(Precision, Precision, Precision)> fn =
            [this, &xy_coords](Precision r, Precision s, Precision /*t*/) -> StaticMatrix<N, N> {
            ShapeFunction   H    = this->shape_function(r, s);      // (N)
            ShapeDerivative dH   = this->shape_derivative(r, s);    // (N×2), dN/dr, dN/ds
            Jacobian        J    = this->jacobian(dH, const_cast<LocalCoords&>(xy_coords));
            Precision       detJ = J.determinant();
            return (H * H.transpose()) * detJ;    // (N×N)
        };
        return this->integration_scheme().integrate(fn);
    }

    MapMatrix mass(Precision* buffer) override {
        // Lokale Achsen & XY-Koordinaten wie bei der Steifigkeit
        auto axes      = get_xyz_axes();
        auto xy_coords = get_xy_coords(axes);

        // Material-/Sektionsdaten
        const Precision rho = this->get_material()->get_density();
        const Precision h   = this->get_section()->thickness;
        const Precision mt  = rho * h;                     // Translationsmasse pro Fläche
        const Precision mr  = rho * (h * h * h / 12.0);    // Rotations-„Masse“ (Flächenträgheit)

        // Flächenintegral der NN^T-Matrix
        StaticMatrix<N, N> M_NN = integrate_NNt(xy_coords);

        // Lokale (x,y,z, rx,ry,rz) konsistente Massenmatrix (6N × 6N)
        StaticMatrix<6 * N, 6 * N> Mloc;
        Mloc.setZero();

        // 1) Translationsblock: ⊗ I3 auf (ux,uy,uz)
        for (Index i = 0; i < N; ++i) {
            for (Index j = 0; j < N; ++j) {
                const Precision mij = mt * M_NN(i, j);
                // (ux,uy,uz) sitzen auf DOFs 0,1,2 im 6er Paket
                for (int k = 0; k < 3; ++k) {
                    const Index I = 6 * i + k;
                    const Index J = 6 * j + k;
                    Mloc(I, J) += mij;
                }
            }
        }

        // 2) Rotationsblock für (rx,ry): ⊗ I2
        for (Index i = 0; i < N; ++i) {
            for (Index j = 0; j < N; ++j) {
                const Precision mij = mr * M_NN(i, j);
                // rx -> dof 3, ry -> dof 4
                for (int k = 0; k < 2; ++k) {
                    const Index I = 6 * i + 3 + k;
                    const Index J = 6 * j + 3 + k;
                    Mloc(I, J) += mij;
                }
            }
        }

        // 3) (Optional) sehr kleiner Drill-Massenanteil für rz (dof 5), damit numerisch stabil
        //    Hier nehmen wir z.B. 1e-6 der mittleren Knotentranslationsmasse.
        {
            Precision avg_m_node = 0.0;
            for (Index i = 0; i < N; ++i)
                avg_m_node += mt * M_NN(i, i);
            avg_m_node /= std::max<Precision>(1.0, static_cast<Precision>(N));
            const Precision m_drill = 1e-6 * avg_m_node;
            for (Index i = 0; i < N; ++i) {
                const Index I = 6 * i + 5;    // rz
                Mloc(I, I) += m_drill;
            }
        }

        // In globale Koordinaten drehen (wie bei der Steifigkeit)
        MapMatrix mapped(buffer, 6 * N, 6 * N);
        auto      T = transformation(axes);    // (6N×6N), setzt pro 3er-Block die Achsen
        mapped      = T.transpose() * Mloc * T;
        return mapped;
    }

	Precision volume() override {
    	// Dicke aus der Section
    	const Precision h = this->get_section()->thickness;

        // Globale Knotentabelle (POSITION) aus dem Modell
        logging::error(this->_model_data != nullptr, "no model data assigned to element ", this->elem_id);
        logging::error(this->_model_data->positions != nullptr, "positions field not set in model data");
        const auto& node_coords_system = *this->_model_data->positions;

    	// Flächeninhalt über die Surface-Geometrie (macht 3D-Jacobian×cross-Produkt)
    	const Precision A = geometry.area(node_coords_system);

    	// Volumen = Dicke * Fläche
    	return h * A;
	}


    virtual Stress stress(Field& displacement, Vec3& xyz) {
        Precision r = xyz(0);
        Precision s = xyz(1);
        Precision t = xyz(2);

        Vec6      res;
        res = Vec6::Zero();

        // get the displacement vectors
        StaticVector<2 * N> disp_membrane;
        StaticVector<3 * N> disp_bending;
        StaticVector<3 * N> disp_shear;

        // get the axes for transformation from global to local
        Mat3 axes = get_xyz_axes();

        for (int i = 0; i < N; i++) {
            ID   node_id             = this->nodes()[i];

            Vec6 displacement_glob   = displacement.row_vec6(static_cast<Index>(node_id));
            Vec3 disp_xyz            = displacement_glob.head(3);
            Vec3 disp_rot            = displacement_glob.tail(3);

            disp_xyz = axes * disp_xyz;   // global → local
            disp_rot = axes * disp_rot;   // global → local

            disp_membrane(2 * i)     = disp_xyz(0);
            disp_membrane(2 * i + 1) = disp_xyz(1);

            disp_bending(3 * i)      = disp_xyz(2);
            disp_bending(3 * i + 1)  = disp_rot(0);
            disp_bending(3 * i + 2)  = disp_rot(1);

            disp_shear(3 * i)        = disp_xyz(2);
            disp_shear(3 * i + 1)    = disp_rot(0);
            disp_shear(3 * i + 2)    = disp_rot(1);
        }

        Mat3 mat_membrane = this->get_elasticity()->get_memb();
        Mat3 mat_bend     = this->get_elasticity()->get_bend(this->get_section()->thickness);
        Mat2 mat_shear    = this->get_elasticity()->get_shear(this->get_section()->thickness);

		// scale material matrices by topo stiffness
        Precision topo_scale = Precision(1);
        if (this->_model_data && this->_model_data->element_stiffness_scale) {
            auto scale_field = this->_model_data->element_stiffness_scale;
            logging::error(scale_field->components == 1,
                           "Field '", scale_field->name, "': element stiffness scale requires 1 component");
            topo_scale = (*scale_field)(static_cast<Index>(this->elem_id));
        }
		mat_membrane *= topo_scale;
		mat_bend     *= topo_scale;
		mat_shear    *= topo_scale;

        // membrane stress
        LocalCoords     xy_coords  = get_xy_coords(axes);
        ShapeFunction   shape_func = this->shape_function(r, s);
        ShapeDerivative shape_der  = this->shape_derivative(r, s);
        Jacobian        jac        = this->jacobian(shape_der, xy_coords);

        // strain displacement matrices
        auto B_membrane = this->strain_disp_membrane(shape_der, jac);
        auto B_bending  = this->strain_disp_bending(shape_der, jac);
        auto B_shear    = this->strain_disp_shear(shape_func, shape_der, jac);

        // membrane stress, s_x, s_y, s_xy
        Vec3 stress_membrane = mat_membrane * B_membrane * disp_membrane;
        res(0) += stress_membrane(0);
        res(1) += stress_membrane(1);
        res(5) += stress_membrane(2);

        // bending stress,
        // mat bend is scaled by h^3 / 12, we need to get rid of that
        Precision h = this->get_section()->thickness;
        mat_bend *= 12 / std::pow(h, 3);
        Precision t_pos          = t * h / 2;
        Vec3      stress_bending = t_pos * mat_bend * B_bending * disp_bending;
        res(0) += stress_bending(0);
        res(1) += stress_bending(1);
        res(5) += stress_bending(2);

        // shear stress
        mat_shear /= this->get_section()->thickness;
        Vec2 stress_shear = mat_shear * B_shear * disp_shear;
        res(3) += stress_shear(0);
        res(4) += stress_shear(1);

        return Stress {res};
    }

    void compute_stress_strain_nodal(Field& displacement,
                                     Field& stress,
                                     Field& strain) override {
        (void) strain; // still unused for now

        // --- Precompute axes and transformation ---
        Mat3 axes   = get_xyz_axes();
        Mat3 axes_T = axes.transpose(); // local -> global

        // --- Precompute thickness + material matrices ---
        Precision h        = this->get_section()->thickness;
        Precision topo_scale = Precision(1);
        if (this->_model_data && this->_model_data->element_stiffness_scale) {
            auto scale_field = this->_model_data->element_stiffness_scale;
            logging::error(scale_field->components == 1,
                           "Field '", scale_field->name, "': element stiffness scale requires 1 component");
            topo_scale = (*scale_field)(static_cast<Index>(this->elem_id));
        }

        Mat3 mat_membrane = this->get_elasticity()->get_memb();
        Mat3 mat_bend     = this->get_elasticity()->get_bend(h);
        Mat2 mat_shear    = this->get_elasticity()->get_shear(h);

        mat_membrane *= topo_scale;
        mat_bend     *= topo_scale;
        mat_shear    *= topo_scale;

        // remove h^3/12 from bending, divide shear by thickness once
        mat_bend  *= 12.0 / (h * h * h);
        mat_shear /= h;

        // --- Build displacement vectors ONCE (global -> local) ---
        StaticVector<2 * N> disp_membrane;
        StaticVector<3 * N> disp_bending;
        StaticVector<3 * N> disp_shear;

        disp_membrane.setZero();
        disp_bending.setZero();
        disp_shear.setZero();

        for (int a = 0; a < N; ++a) {
            ID   node_id           = this->nodes()[a];

            Vec6 displacement_glob = displacement.row_vec6(static_cast<Index>(node_id));
            Vec3 disp_xyz          = displacement_glob.head(3);
            Vec3 disp_rot          = displacement_glob.tail(3);

            // global -> local only once
            disp_xyz = axes * disp_xyz;
            disp_rot = axes * disp_rot;

            disp_membrane(2 * a    ) = disp_xyz(0);
            disp_membrane(2 * a + 1) = disp_xyz(1);

            disp_bending(3 * a    ) = disp_xyz(2);
            disp_bending(3 * a + 1) = disp_rot(0);
            disp_bending(3 * a + 2) = disp_rot(1);

            // your current MITC-like shear uses same dofs as bending
            disp_shear(3 * a    ) = disp_xyz(2);
            disp_shear(3 * a + 1) = disp_rot(0);
            disp_shear(3 * a + 2) = disp_rot(1);
        }

        // --- Geometry: local node coordinates & xy coords once ---
        StaticMatrix<N, 2> coords    = this->geometry.node_coords_local();
        LocalCoords        xy_coords = get_xy_coords(axes);

        // --- Loop over nodes and evaluate stress at top/bottom ---
        for (int i = 0; i < N; ++i) {
            ID node_id = this->nodes()[i];

            Precision r = coords(i, 0);
            Precision s = coords(i, 1);

            // shape & jacobian at this (r,s): same for top/bottom
            ShapeFunction   shape_func = this->shape_function(r, s);
            ShapeDerivative shape_der  = this->shape_derivative(r, s);
            Jacobian        jac        = this->jacobian(shape_der, xy_coords);

            auto B_membrane = this->strain_disp_membrane(shape_der, jac);
            auto B_bending  = this->strain_disp_bending(shape_der, jac);
            auto B_shear    = this->strain_disp_shear(shape_func, shape_der, jac);

            // mid-surface membrane + shear stress: independent of t
            Vec6 res_mid = Vec6::Zero();

            Vec3 stress_membrane = mat_membrane * (B_membrane * disp_membrane);
            res_mid(0) += stress_membrane(0);
            res_mid(1) += stress_membrane(1);
            res_mid(5) += stress_membrane(2);

            Vec2 stress_shear = mat_shear * (B_shear * disp_shear);
            res_mid(3) += stress_shear(0);
            res_mid(4) += stress_shear(1);

            // bending contribution: changes sign with t
            Vec6 res_bend_top = Vec6::Zero();
            Vec6 res_bend_bot = Vec6::Zero();

            Vec3 kappa      = B_bending * disp_bending; // curvature vector
            Vec3 sigma_bend = mat_bend * kappa;         // bending stress factor

            Precision t_top = +1.0;
            Precision t_bot = -1.0;
            Precision z_top = t_top * h * 0.5;
            Precision z_bot = t_bot * h * 0.5;

            // add bending part
            res_bend_top(0) += z_top * sigma_bend(0);
            res_bend_top(1) += z_top * sigma_bend(1);
            res_bend_top(5) += z_top * sigma_bend(2);

            res_bend_bot(0) += z_bot * sigma_bend(0);
            res_bend_bot(1) += z_bot * sigma_bend(1);
            res_bend_bot(5) += z_bot * sigma_bend(2);

            Stress stress_top_local{ res_mid + res_bend_top };
            Stress stress_bot_local{ res_mid + res_bend_bot };

            // pick the one with larger norm (still in local coords)
            Stress stress_nodal_local = (stress_top_local.norm() > stress_bot_local.norm())
                                        ? stress_top_local
                                        : stress_bot_local;

            // transform to global once with axes_T
            Stress stress_nodal_global = stress_nodal_local.transform(axes_T);

            const Index node_idx = static_cast<Index>(node_id);
            for (int j = 0; j < 6; ++j) {
                stress(node_idx, j) += stress_nodal_global(j);
            }
        }
    }


    // in DefaultShellElement<...>
    void compute_stress_strain(Field& ip_stress,
                               Field& /*ip_strain*/,
                               Field& displacement,
                               int       ip_offset) override {
        Mat3            axes      = get_xyz_axes();
        auto            xy_coords = get_xy_coords(axes);
        const Precision h         = this->get_section()->thickness;

        // Membrane material (σx, σy, τxy)
        Mat3 Dm = this->get_elasticity()->get_memb();

        // Local membrane displacement vector [ux, uy] per node
        StaticVector<2 * N> u_mem;
        u_mem.setZero();
        for (Index i = 0; i < N; ++i) {
            ID   node_id     = this->nodes()[i];
            Vec6 u_glob      = displacement.row_vec6(static_cast<Index>(node_id));    // (ux,uy,uz,rx,ry,rz)
            Vec3 t_glob      = u_glob.head<3>();
            Vec3 t_loc       = axes * t_glob;
            u_mem(2 * i)     = t_loc(0);    // ux
            u_mem(2 * i + 1) = t_loc(1);    // uy
        }

        const auto& scheme = this->integration_scheme();
        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const Precision r     = scheme.get_point(ip).r;
            const Precision s     = scheme.get_point(ip).s;

            ShapeDerivative dH_rs = this->shape_derivative(r, s);
            Jacobian        J     = this->jacobian(dH_rs, xy_coords);

            // B-matrix for membrane part (3×2N)
            StaticMatrix<3, 2 * N> Bm = this->strain_disp_membrane(dH_rs, J);

            // Membrane stresses in local coords
            StaticVector<3> sigma = Dm * (Bm * u_mem);

            // Stress resultants (force per unit length)
            const Precision nxx = sigma(0) * h;
            const Precision nyy = sigma(1) * h;
            const Precision nxy = sigma(2) * h;

            const Index     row = static_cast<Index>(ip_offset) + ip;
            ip_stress(row, 0)   = nxx;    // xx
            ip_stress(row, 1)   = nyy;    // yy
            ip_stress(row, 2)   = 0.0;    // zz
            ip_stress(row, 3)   = 0.0;    // yz
            ip_stress(row, 4)   = 0.0;    // zx
            ip_stress(row, 5)   = nxy;    // xy
        }
    }
};

}    // namespace fem::model

#endif    // SHELL_SIMPLE_H
