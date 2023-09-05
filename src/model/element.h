#pragma once

#include "../core/core.h"
#include "../material/material.h"
#include "quadrature.h"

#include <array>

namespace fem {

namespace model{

struct ElementInterface{
    const ID            elem_id  = 0;
    protected:
    material::Material* material = nullptr;

    public:
    ElementInterface(ID p_elem_id)
        : elem_id(p_elem_id) {}

    void set_material(material::Material* material) {
        ElementInterface::material = material;
    }

    virtual MapMatrix stiffness(Precision* buffer) = 0;
    virtual Dofs dofs() = 0;
    virtual Dim dimensions() = 0;
    virtual Dim n_nodes() = 0;
    virtual ID* nodes() = 0;

    virtual void compute_stress    (NodeData& node_coords, NodeData& displacement, NodeData& stress, NodeData& strain) = 0;
    virtual void compute_compliance(NodeData& node_coords, NodeData& displacement, ElementData& result) = 0;
};

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
    virtual StaticMatrix<n_strain, D * N> strain_displacement (const StaticMatrix<N, D> &shape_der_global) = 0;
    virtual StaticMatrix<n_strain, D * N> strain_displacements(const StaticMatrix<N, D> &node_coords, Precision r, Precision s, Precision t, Precision& det){
        StaticMatrix<N, D> local_shape_der = shape_derivative(r, s, t);
        StaticMatrix<D, D> jac             = jacobian(node_coords, r, s, t);

        det = jac.determinant();
        StaticMatrix<D, D> inv = jac.inverse();

        logging::error(det > 0,
                  "negative determinant encountered in element ",elem_id," det: ",det);
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
    MapMatrix stiffness(Precision* buffer) override {
        logging::error(material != nullptr, "no material assigned to element ", elem_id);
        logging::error(material->has_elasticity(), "material has no elasticity components assigned at element ", elem_id);

        StaticMatrix<N, D> node_coords = this->node_coords_local();

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


    void compute_stress    (NodeData& node_coords, NodeData& displacement, NodeData& stress, NodeData& strain) override {
        auto local_node_coords  = this->node_coords_local();

        auto local_disp_mat     = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
        auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

        for(int n = 0; n < n_nodes(); n++){

            Precision r = local_node_coords(n, 0);
            Precision s = local_node_coords(n, 0);
            Precision t = local_node_coords(n, 0);
            Precision det;

            StaticMatrix<n_strain, D * N>    B = this->strain_displacements(local_node_coords, r, s, t, det);
            StaticMatrix<n_strain, n_strain> E = this->material->elasticity()->template get<D>();

            auto node_id = node_ids[n];

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
        auto K = stiffness(buffer);

        // 2. Extract and flatten nodal displacements
        auto local_disp_mat     = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
        auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

        Precision strain_energy = local_displacement.dot((K * local_displacement));

        if(strain_energy < -0.01){
            Eigen::EigenSolver< DynamicMatrix > solver{DynamicMatrix(K)};
            std::cout << "====================" << std::endl;
            std::cout << solver.eigenvalues() << std::endl;
            std::cout << (local_displacement).transpose() << std::endl;
            std::cout << (K * local_displacement).transpose() << std::endl;
            std::cout << (K)  << std::endl;
            std::cout << strain_energy  << std::endl;
            std::cout << "====================" << std::endl;
            exit(0);
        }

        result(elem_id, 0) = strain_energy;
    }

};

using ElementPtr = std::shared_ptr<ElementInterface>;

#include "element.ipp"

}

}    // namespace fem
