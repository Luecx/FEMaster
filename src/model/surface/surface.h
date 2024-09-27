#pragma once

#include "../element.h"
#include "../element_solid.h"
#include "../../math/quadrature.h"
#include <functional>

namespace fem::model {

/******************************************************************************
 * @class SurfaceInterface
 * @brief Base class template for a surface element in finite element analysis.
 *
 * @tparam N The number of nodes in the surface element.
 ******************************************************************************/
template<Index N>
struct SurfaceInterface {

    std::array<ID, N> nodeIds;

    SurfaceInterface(const std::array<ID, N>& pNodeIds) : nodeIds(pNodeIds) {}

    // ------------
    virtual StaticMatrix<N, 1> shape_function  (Precision r, Precision s) const = 0;
    virtual StaticMatrix<N, 2> shape_derivative(Precision r, Precision s) const = 0;

    virtual StaticMatrix<3, 2> jacobian(const StaticMatrix<N, 3>& node_coords, Precision r, Precision s) const {
        StaticMatrix<N, 2> local_shape_derivative = shape_derivative(r, s);
        StaticMatrix<3, 2> jacobian {};

        // Compute the 3x2 Jacobian matrix, which maps local (r, s) to 3D global space
        for (Dim m = 0; m < 2; m++) { // Derivatives with respect to r, s
            for (Dim n = 0; n < 3; n++) { // x, y, z global coordinates
                Precision dxn_drm = 0;
                for (Dim k = 0; k < N; k++) { // Sum over all nodes
                    dxn_drm += node_coords(k, n) * local_shape_derivative(k, m);
                }
                jacobian(n, m) = dxn_drm;
            }
        }

        return jacobian;
    }

    virtual StaticMatrix<N, 2> node_coords_local() const = 0;

    virtual StaticMatrix<N, 3> node_coords_global(const NodeData& node_coords_system) const {
        StaticMatrix<N, 3> res {};
        for (Index i = 0; i < N; i++) {
            for(Index j = 0; j < 3; j++) {
                res(i, j) = node_coords_system(nodeIds[i], j);
            }
        }
        return res;
    }

    // ------------
    virtual StaticVector<3> local_to_global(const StaticVector<2>& local, const NodeData& node_coords_system) const {
        auto node_coords_global = this->node_coords_global(node_coords_system);
        StaticVector<3> res{};
        for (Index i = 0; i < N; i++) {
            res += node_coords_global.row(i) * shape_function(local(0), local(1))(i);
        }
        return res;
    }

    // projection
    virtual StaticVector<2> global_to_local(const StaticVector<3>& global, const NodeData& node_coords_system, bool clip) const = 0;

    // ------------
    virtual const quadrature::Quadrature& integration_scheme() const = 0;

    Precision min_node_distance(const StaticVector<3>& nodeCoord, const NodeData& node_coords_system) const {
        // compute the area distance
        Precision planar_distance = (local_to_global(global_to_local(nodeCoord, node_coords_system), node_coords_system) - nodeCoord).norm();
        return planar_distance;
    }

    Precision integrate_constant(Precision value, const NodeData& node_coords_system) const {
        return this->area(node_coords_system) * value;
    }

    Precision area(const NodeData& node_coords_system) const {
        auto node_coords_global = this->node_coords_global(node_coords_system);

        // Define the integrand for computing surface area
        std::function<Precision(Precision, Precision, Precision)> integrand = [&](Precision r, Precision s, Precision) -> Precision {
            // Compute the 3x2 Jacobian
            auto jac = jacobian(node_coords_global, r, s);

            // Compute the cross product of the Jacobian columns
            StaticVector<3> tangent_r = jac.col(0); // ∂r/∂ξ
            StaticVector<3> tangent_s = jac.col(1); // ∂r/∂η

            // Compute the magnitude of the cross product: |tangent_r x tangent_s|
            Precision area_element = (tangent_r.cross(tangent_s)).norm();

            return area_element;
        };

        // Perform numerical integration using the specified integration scheme
        Precision total_area = integration_scheme().integrate(integrand);
        return total_area;
    }

    const std::array<ID, N>& get_node_ids() const {
        return nodeIds;
    }
};

} // namespace fem::model
