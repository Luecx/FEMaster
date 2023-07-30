#pragma once

#include "../core/core.h"
#include "../material/material.h"
#include "quadrature.h"

#include <array>

namespace fem {

namespace model{

struct ElementInterface{
    protected:
    const ID              elem_id  = 0;
    material::MaterialPtr material = nullptr;

    public:
    ElementInterface(ID p_elem_id)
        : elem_id(p_elem_id) {}

    void set_material(const material::MaterialPtr& material) {
        ElementInterface::material = material;
    }

    virtual DynamicMatrix stiffness() = 0;
    virtual Dofs dofs() = 0;
    virtual Dim dimensions() = 0;
    virtual Dim n_nodes() = 0;
    virtual ID* nodes() = 0;
};

template<size_t N, Dim D>
struct SolidElement : public ElementInterface{
    public:
    std::array<ID, N> node_ids{};

    protected:
    static constexpr Dim n_strain = D == 2 ? 3:6;

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
    virtual StaticMatrix<n_strain, D * N> strain_displacement(const StaticMatrix<N, D> &shape_der_global) = 0;
    virtual StaticMatrix<n_strain, D * N> strain_displacements(const StaticMatrix<N, D> &node_coords, Precision r, Precision s, Precision t, Precision& det){
        StaticMatrix<N, D> local_shape_der = shape_derivative(r, s, t);
        StaticMatrix<D, D> jac             = jacobian(node_coords, r, s, t);

        det = jac.determinant();
        StaticMatrix<D, D> inv = jac.inverse();

        log_error(det > 0,
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
        return Dofs {{true}, {true}, {D > 2}, {false}, {false}, {false}};
    }
    DynamicMatrix stiffness() override {
        log_error(material != nullptr, "no material assigned to element ", elem_id);
        log_error(material->has_elasticity(), "material has no elasticity components assigned at element ", elem_id);

        StaticMatrix<N, D> node_coords = this->node_coords_local();

        std::function<StaticMatrix<D * N, D * N>(Precision, Precision, Precision)> func = [this, node_coords](Precision r, Precision s, Precision t) {
                Precision                        det;
                StaticMatrix<n_strain, D* N>     B   = this->strain_displacements(node_coords, r, s, t, det);
                StaticMatrix<n_strain, n_strain> E   = this->material->elasticity()->template get<D>();
                StaticMatrix<D * N, D* N>        res = B.transpose() * (E * B) * det;
                return StaticMatrix<D * N, D* N>(res);
            };
        StaticMatrix<D * N, D * N> stiffness = integration_scheme().integrate(func);
        return DynamicMatrix(stiffness);
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
};

using ElementPtr = std::shared_ptr<ElementInterface>;

#include "element.ipp"

}

}    // namespace fem
