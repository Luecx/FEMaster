#pragma once

#include "element.h"

namespace fem::model {

template<size_t N>
struct SolidElement : public ElementInterface{
    public:
    std::array<ID, N> node_ids{};

    protected:
    static constexpr Dim n_strain = 6;
    static constexpr Dim D = 3;

    friend quadrature::Quadrature;

    public:
    SolidElement(ID p_elem_id, std::array<ID, N> p_node_ids)
        : ElementInterface(p_elem_id)
        , node_ids(p_node_ids) {}

    // shape function and its derivative, note that the derivative is in local coordinates
    // for transforming into global coordinates, the inverse of the jacobian is required
    virtual StaticMatrix<N, 1> shape_function  (Precision r, Precision s, Precision t) = 0;
    virtual StaticMatrix<N, D> shape_derivative(Precision r, Precision s, Precision t) = 0;

    // node coords in local coordinates and global coordinates
    virtual StaticMatrix<N, D> node_coords_local() = 0;
    virtual StaticMatrix<N, D> node_coords_global(NodeData& node_coords) {
        return this->nodal_data<D>(node_coords);
    }

    // integration scheme
    virtual const quadrature::Quadrature& integration_scheme() = 0;

    // compute the strain displacement matrix B from the compact form
    virtual StaticMatrix<n_strain, D * N> strain_displacement (const StaticMatrix<N, D> &shape_der_global) {
        StaticMatrix<n_strain, D * N> B {};
        B.setZero();
        for (int j = 0; j < N; j++) {
            int r1   = j  * 3;
            int r2   = r1 + 1;
            int r3   = r1 + 2;
            B(0, r1) = shape_der_global(j, 0);
            B(1, r2) = shape_der_global(j, 1);
            B(2, r3) = shape_der_global(j, 2);
            B(3, r1) = shape_der_global(j, 1) / (Precision) 2.0;    // divide by 2 in order to account for real shear strain
            B(3, r2) = shape_der_global(j, 0) / (Precision) 2.0;    // divide by 2 in order to account for real shear strain
            B(4, r1) = shape_der_global(j, 2) / (Precision) 2.0;    // divide by 2 in order to account for real shear strain
            B(4, r3) = shape_der_global(j, 0) / (Precision) 2.0;    // divide by 2 in order to account for real shear strain
            B(5, r2) = shape_der_global(j, 2) / (Precision) 2.0;    // divide by 2 in order to account for real shear strain
            B(5, r3) = shape_der_global(j, 1) / (Precision) 2.0;    // divide by 2 in order to account for real shear strain
        }
        return B;
    }
    virtual StaticMatrix<n_strain, D * N> strain_displacements(const StaticMatrix<N, D> &node_coords, Precision r, Precision s, Precision t, Precision& det){
        StaticMatrix<N, D> local_shape_der = shape_derivative(r, s, t);
        StaticMatrix<D, D> jac             = jacobian(node_coords, r, s, t);

        det = jac.determinant();
        StaticMatrix<D, D> inv = jac.inverse();
        logging::error(det > 0,"negative determinant encountered in element ",elem_id,
                                "\ndet        : ", det,
                                "\nCoordinates: ", node_coords,
                                "\nJacobi     : ", jac);
        StaticMatrix<N, D> global_shape_der = (inv * local_shape_der.transpose()).transpose();

        return strain_displacement(global_shape_der);
    }

    // get the jacobian -> derivatives of x,y,z with respect to r,s,t
    // requires the coordinates of the nodes. this can be retrieved using nodal_data from the global node location
    // matrix
    virtual StaticMatrix<D, D> jacobian(const StaticMatrix<N, D>& node_coords, Precision r, Precision s, Precision t){
        StaticMatrix<N, D> local_shape_derivative = shape_derivative(r, s, t);
        StaticMatrix<D, D> jacobian {};

        for(int m = 0; m < D; m++){
            for(int n = 0; n < D; n++){
                Precision dxn_drm = 0;
                for(int k = 0; k < N; k++){
                    dxn_drm += node_coords(k, n) * local_shape_derivative(k,m);
                }
                jacobian(m, n) = dxn_drm;
            }
        }

        return jacobian;
    }

    template<size_t K>
    StaticMatrix<N, K> nodal_data(const NodeData &full_data, int offset=0, int stride=1){
        StaticMatrix<N, K> res {};
        runtime_assert(full_data.cols() >= offset + stride * D, "cannot extract this many elements from the data");
        for (int m = 0; m < N; m++) {
            for (int j = 0; j < D; j++) {
                int n     = j * stride + offset;
                res(m, j) = full_data(node_ids[m], n);
            }
        }
        return res;
    }

    public:
    Dofs dofs() override {
        // degrees of freedom for each node is dx,dy,dz
        // no dof for rx,ry,rz
        return Dofs {{true}, {true}, {true}, {false}, {false}, {false}};
    }
    MapMatrix stiffness(NodeData& position, Precision* buffer) override {
        logging::error(material != nullptr, "no material assigned to element ", elem_id);
        logging::error(material->has_elasticity(), "material has no elasticity components assigned at element ", elem_id);

        StaticMatrix<N, D> node_coords = this->node_coords_global(position);

        std::function<StaticMatrix<D * N, D * N>(Precision, Precision, Precision)> func =
            [this, node_coords](Precision r, Precision s, Precision t) {
                Precision                        det;
                StaticMatrix<n_strain, D * N>    B   = this->strain_displacements(node_coords, r, s, t, det);
                StaticMatrix<n_strain, n_strain> E   = this->material->elasticity()->template get<D>();
                StaticMatrix<D * N, D * N>        res = B.transpose() * (E * B) * det;
                return StaticMatrix<D * N, D* N>(res);
            };
        StaticMatrix<D * N, D * N> stiffness = integration_scheme().integrate(func);
        stiffness = 0.5 * (stiffness + stiffness.transpose());

        MapMatrix mapped{buffer, D * N, D * N};
        mapped = stiffness;
        return mapped;
    }

    Dim dimensions() override {
        return D;
    }

    Dim n_nodes() override {
        return node_ids.size();
    }

    ID* nodes() override {
        return &node_ids[0];
    }

    Precision volume(NodeData& node_coords) override {
        StaticMatrix<N, D> node_coords_glob = this->node_coords_global(node_coords);
        std::function<Precision(Precision, Precision, Precision)> func =
            [this, node_coords_glob](Precision r, Precision s, Precision t) {
                Precision det = jacobian(node_coords_glob, r, s, t).determinant();
                return det;
            };
        Precision volume = integration_scheme().integrate(func);
        return volume;
    }

    void compute_stress    (NodeData& node_coords, NodeData& displacement, NodeData& stress, NodeData& strain) override {
        auto local_node_coords  = this->node_coords_local();
        auto global_node_coords = this->node_coords_global(node_coords);

        auto local_disp_mat     = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
        auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

        for(int n = 0; n < n_nodes(); n++){

            Precision r = local_node_coords(n, 0);
            Precision s = local_node_coords(n, 1);
            Precision t = local_node_coords(n, 2);
            Precision det;

            StaticMatrix<n_strain, D * N>    B = this->strain_displacements(global_node_coords, r, s, t, det);
            StaticMatrix<n_strain, n_strain> E = this->material->elasticity()->template get<D>();

            auto node_id  = node_ids[n];
            auto strains  = B * local_displacement;
            auto stresses = E * strains;

            for(int j = 0; j < n_strain; j++){
                strain(node_id, j) += strains(j);
                stress(node_id, j) += stresses(j);
            }

        }
    }
    void compute_compliance(NodeData& node_coords, NodeData& displacement, ElementData& result) override {
        // 1. Compute the element stiffness matrix
        Precision buffer[D * N * D * N];
        auto K = stiffness(node_coords, buffer);

        // 2. Extract and flatten nodal displacements
        auto local_disp_mat     = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
        auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

        Precision strain_energy = local_displacement.dot((K * local_displacement));
        result(elem_id, 0) = strain_energy;
    }

    template<class ElementType>
    static void test_implementation() {
        // Create an instance of the element type
        std::array<ID, N> nodeArray;
        for (size_t i = 0; i < N; i++) {
            nodeArray[i] = static_cast<ID>(i);    // Just initializing to example node IDs
        }

        ElementType el(0, nodeArray);

        // Test 1: Checking node values
        auto               node_coords = el.node_coords_local();
        StaticMatrix<N, N> globalMatrix;
        globalMatrix.setZero();

        for (size_t i = 0; i < N; i++) {
            Precision          r               = node_coords(i, 0);
            Precision          s               = node_coords(i, 1);
            Precision          t               = node_coords(i, 2);

            StaticMatrix<N, 1> shapeFuncValues = el.shape_function(r, s, t);
            for (size_t j = 0; j < N; j++) {
                globalMatrix(i, j) = shapeFuncValues(j);
            }
        }
        std::cout << globalMatrix << std::endl;

        // Test 2: Checking shape function sum
        const Precision step      = 0.2;
        const Precision tolerance = 1e-6;

        for (Precision r = -1; r <= 1; r += step) {
            for (Precision s = -1; s <= 1; s += step) {
                for (Precision t = -1; t <= 1; t += step) {
                    StaticMatrix<N, 1> shapeFuncValues = el.shape_function(r, s, t);

                    Precision          sum             = 0;
                    for (size_t j = 0; j < N; j++) {
                        sum += shapeFuncValues(j);
                    }

                    if (std::abs(sum - 1.0) > tolerance) {
                        std::cout << "Sum of shape functions at (r, s, t) = (" << r << ", " << s << ", " << t
                                  << ") is not 1. Actual sum: " << sum << std::endl;
                    }
                }
            }
        }

        const Precision delta = 1e-6;
        for (Precision r = -1; r <= 1; r += step) {
            for (Precision s = -1; s <= 1; s += step) {
                for (Precision t = -1; t <= 1; t += step) {

                    // Derivative from the function
                    StaticMatrix<N, D> true_derivatives = el.shape_derivative(r, s, t);

                    // Compute finite differences for each direction
                    StaticMatrix<N, 1> shapeFuncValues_r_plus_delta = el.shape_function(r + delta, s, t);
                    StaticMatrix<N, 1> shapeFuncValues_s_plus_delta = el.shape_function(r, s + delta, t);
                    StaticMatrix<N, 1> shapeFuncValues_t_plus_delta = el.shape_function(r, s, t + delta);

                    StaticMatrix<N, 1> shapeFuncValues_r_minu_delta = el.shape_function(r - delta, s, t);
                    StaticMatrix<N, 1> shapeFuncValues_s_minu_delta = el.shape_function(r, s - delta, t);
                    StaticMatrix<N, 1> shapeFuncValues_t_minu_delta = el.shape_function(r, s, t - delta);

                    StaticMatrix<N, 1> shapeFuncValues = el.shape_function(r, s, t);

                    StaticMatrix<N, D> finite_diff_derivatives;

                    for (size_t j = 0; j < N; j++) {
                        finite_diff_derivatives(j, 0) = (shapeFuncValues_r_plus_delta(j) - shapeFuncValues_r_minu_delta(j)) / (2 * delta);  // dr
                        finite_diff_derivatives(j, 1) = (shapeFuncValues_s_plus_delta(j) - shapeFuncValues_s_minu_delta(j)) / (2 * delta);  // ds
                        finite_diff_derivatives(j, 2) = (shapeFuncValues_t_plus_delta(j) - shapeFuncValues_t_minu_delta(j)) / (2 * delta);  // dt
                    }

                    // Compare true derivatives with finite differences
                    for (size_t j = 0; j < N; j++) {
                        for (size_t d = 0; d < D; d++) {
                            if (std::abs(true_derivatives(j, d) - finite_diff_derivatives(j, d)) > tolerance) {
                                std::cout << "Mismatch in derivative at (r, s, t) = (" << r << ", " << s << ", " << t
                                          << ") in direction " << d
                                          << " from shape function " << j
                                          << ". True derivative: " << true_derivatives(j, d)
                                          << ", Finite Difference: " << finite_diff_derivatives(j, d) << std::endl;
                            }
                        }
                    }
                }
            }
        }

    }

};


}