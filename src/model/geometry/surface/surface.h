/**
 * @file surface.hpp
 * @brief Defines the SurfaceInterface template class, which serves as a base
 *        class for surface finite elements in both 2D and 3D spaces. It provides
 *        methods for calculating shape functions, projecting points onto the surface,
 *        computing local coordinates, and integrating functions across the surface.
 *
 * @tparam N The number of nodes associated with the surface element.
 *
 * @details This class provides a generic interface for surface elements in FEM.
 *          It defines functions for computing shape functions and derivatives,
 *          projecting points in 3D space to the surface, computing local coordinates
 *          for projected points, and performing surface integration using a lambda
 *          function. Derived classes should implement the specific shape functions
 *          for different surface element types (e.g., triangular, quadrilateral,
 *          first- or second-order).
 *
 *          The surface types supported include triangular and quadrilateral elements
 *          with varying order (e.g., S3D3, S3D4, S3D6, S3D8), and the interface provides
 *          a framework for computing element-surface interactions and enabling complex
 *          tie constraints in FEM. This is not to be confused with SHELL elements,
 *          as it focuses on surface interactions in a 3D domain.
 *
 * @date Created on 27.09.2024
 * @author Finn Eggers
 */

#pragma once

#include "../../../core/core.h"
#include "../../../math/quadrature.h"

#include <Eigen/Dense>
#include <functional>
#include <iostream>
#include <memory>

namespace fem::model {

struct SurfaceInterface {
    const Index n_edges;
    const Index n_nodes;
    const Index n_nodes_per_edge;

    /**
     * @brief Constructor for SurfaceInterface
     *
     * @param nEdges Number of edges in the surface.
     * @param nNodes Number of nodes in the surface.
     * @param nNodesPerEdge Number of nodes per edge.
     */
    SurfaceInterface(Index nEdges = 0, Index nNodes = 0, Index nNodesPerEdge = 0)
        : n_edges(nEdges)
        , n_nodes(nNodes)
        , n_nodes_per_edge(nNodesPerEdge) {}

    /**
     * @brief Virtual destructor to ensure proper polymorphic behavior.
     */
    virtual ~SurfaceInterface()                                                               = default;

    virtual Vec3 local_to_global(const Vec2& local, const NodeData& node_coords_system) const = 0;
    virtual Vec2 global_to_local(const Vec3& global, const NodeData& node_coords_system, bool clip = false) const = 0;
    virtual Vec3 normal(const NodeData& node_coords_system, const Vec2& local) const                              = 0;
    virtual bool in_bounds(const Vec2& local) const                                                               = 0;
    virtual Precision     area(const NodeData& node_coords_system) const                                          = 0;
    virtual DynamicVector shape_function(const Vec2& local) const                                                 = 0;
    virtual DynamicVector shape_function_integral(const NodeData& node_coords_system) const                       = 0;
    virtual ID*           nodes()                                                                                 = 0;

    virtual void          apply_dload(NodeData& node_coords, NodeData& node_loads, Vec3 load)                     = 0;
    virtual void          apply_pload(NodeData& node_coords, NodeData& node_loads, Precision load)                = 0;

    ID*                   begin() {
        return nodes();
    }
    ID* end() {
        return nodes() + n_nodes;
    }
    using Ptr = std::shared_ptr<SurfaceInterface>;
};

/**
 * @class Surface
 * @brief Base class template for a surface element in finite element analysis.
 *
 * @tparam N The number of nodes in the surface element.
 *
 * @details SurfaceInterface serves as the foundational interface for implementing
 *          surface elements such as triangular and quadrilateral elements. The class
 *          provides virtual methods for calculating shape functions, computing derivatives,
 *          performing projection of 3D points to the surface, and conducting surface integration.
 */
template<Index N>
struct Surface : public SurfaceInterface {
    static constexpr Index num_edges          = (N > 4 ? N / 2 : N);
    static constexpr Index num_nodes          = N;
    static constexpr Index num_nodes_per_edge = N > 4 ? 3 : 2;

    /**
     * @brief Array of node IDs corresponding to the surface element.
     */
    std::array<ID, N> nodeIds;

    /**
     * @brief Constructor for Surface.
     *
     * @param pNodeIds Array of node IDs for the surface element.
     */
    Surface(const std::array<ID, N>& pNodeIds)
        : SurfaceInterface(num_edges, num_nodes, num_nodes_per_edge)
        , nodeIds(pNodeIds) {}

    virtual ~Surface() = default;

    /**
     * @brief Compute the shape functions at a given local coordinate (r, s).
     *
     * @param r Local coordinate in the parametric space of the surface element.
     * @param s Local coordinate in the parametric space of the surface element.
     * @return StaticMatrix<N, 1> Vector of shape function values at (r, s).
     */
    virtual StaticMatrix<N, 1> shape_function(Precision r, Precision s) const = 0;

    /**
     * @brief interpolates a matrix of shape functions at a given local coordinate (r, s).
     **/
    template<Index M>
    StaticVector<M> interpolate(const StaticMatrix<N, M>& nodal_values, Precision r, Precision s) const {
        StaticVector<M> res = (shape_function(r,s).transpose() * nodal_values).transpose();
    }

    /**
     * @brief Gives a dynamic vector for outside use which is a bit slower.
     */
    DynamicVector shape_function(const Vec2& local) const override {
        return DynamicVector(shape_function(local(0), local(1)));
    }

    /**
     * @brief Compute the first derivative of shape functions at a given local coordinate (r, s).
     *
     * @param r Local coordinate in the parametric space of the surface element.
     * @param s Local coordinate in the parametric space of the surface element.
     * @return StaticMatrix<N, 2> Matrix of shape function derivatives at (r, s).
     *
     * @details The returned matrix contains the derivatives with respect to r and s for each node.
     *          Rows correspond to nodes, and columns correspond to derivatives with respect to r and s.
     */
    virtual StaticMatrix<N, 2> shape_derivative(Precision r, Precision s) const = 0;

    /**
     * @brief Compute the second derivative of shape functions at a given local coordinate (r, s).
     *
     * @param r Local coordinate in the parametric space of the surface element.
     * @param s Local coordinate in the parametric space of the surface element.
     * @return StaticMatrix<N, 3> Matrix of second-order shape function derivatives at (r, s).
     *
     * @details The returned matrix contains the second derivatives with respect to r², s², and r*s
     *          for each node. Rows correspond to nodes, and columns correspond to the respective
     *          second-order derivative.
     */
    virtual StaticMatrix<N, 3> shape_second_derivative(Precision r, Precision s) const = 0;

    /**
     * @brief Compute the Jacobian matrix at a given local coordinate (r, s).
     *
     * @param node_coords Node coordinates in the global space.
     * @param r Local coordinate in the parametric space of the surface element.
     * @param s Local coordinate in the parametric space of the surface element.
     * @return StaticMatrix<3, 2> 3x2 Jacobian matrix mapping local coordinates to global coordinates.
     *
     * @details The Jacobian matrix is used for transformations between local and global coordinates.
     *          It is calculated by multiplying the shape function derivatives with the node coordinates.
     */
    StaticMatrix<3, 2> jacobian(const StaticMatrix<N, 3>& node_coords, Precision r, Precision s) const {
        StaticMatrix<N, 2> local_shape_derivative = shape_derivative(r, s);
        StaticMatrix<3, 2> jacobian {};

        // Compute the 3x2 Jacobian matrix, which maps local (r, s) to 3D global space
        for (Dim m = 0; m < 2; m++) {        // Derivatives with respect to r, s
            for (Dim n = 0; n < 3; n++) {    // x, y, z global coordinates
                Precision dxn_drm = 0;
                for (Dim k = 0; k < N; k++) {    // Sum over all nodes
                    dxn_drm += node_coords(k, n) * local_shape_derivative(k, m);
                }
                jacobian(n, m) = dxn_drm;
            }
        }

        return jacobian;
    }

    /**
     * @brief Sets a local coordinate to its corresponding global position on the surface.
     *
     * @param local Vector containing the local coordinates (r, s).
     * @param node_coords_system Node coordinates in the global system.
     * @return Vec3 Global position corresponding to the local coordinate.
     *
     * @details Computes the global coordinates corresponding to a local (r, s) coordinate
     *          by evaluating the shape functions at (r, s) and summing the contributions
     *          from each node's global coordinates.
     */
    Vec3 local_to_global(const Vec2& local, const NodeData& node_coords_system) const override {
        auto node_coords_global = this->node_coords_global(node_coords_system);
        Vec3 res                = Vec3::Zero();
        for (Index i = 0; i < N; i++) {
            res += node_coords_global.row(i) * shape_function(local(0), local(1))(i);
        }
        return res;
    }

    /**
     * @brief Find the local coordinates (r, s) corresponding to the closest point on the surface element to a given
     * global point.
     *
     * @param p The point in global coordinates to which the closest point on the surface element is sought.
     * @param node_coords_system Node coordinates in the global coordinate system.
     * @param clip If true, the function will also consider boundary solutions using line elements when the solution is
     * out of bounds.
     * @return Vec2 The local coordinates (r, s) corresponding to the closest point on the surface element to point 'p'.
     *
     * @details This function uses an iterative method (Newton-Raphson) to find the parametric coordinates (r, s) that
     *          minimize the distance between the point 'p' and the position vector of the surface element. Multiple
     *          initial guesses are used to ensure convergence to the global minimum distance. If the solution is out of
     * bounds and clip is true, the function will check the boundaries using the corresponding line elements.
     */
    Vec2 global_to_local(const Vec3& p, const NodeData& node_coords_system, bool clip = false) const override {
        constexpr int       max_iter = 32;
        constexpr Precision eps      = 1e-12;

        std::vector<Vec2>   initial_guesses;

        // Determine initial guesses based on N
        switch (N) {
            case 3:    // Triangle element with 3 nodes
                initial_guesses = {{0.25, 0.25}};
                break;
            case 4:    // Quadrilateral element with 4 nodes
                initial_guesses = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
                break;
            case 6:    // Triangle element with 6 nodes (quadratic)
                initial_guesses = {{0, 0}, {1, 0}, {0, 1}};
                break;
            case 8:    // Quadrilateral element with 8 nodes (quadratic)
                initial_guesses = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {0, -1}, {1, 0}, {0, 1}, {-1, 0}};
                break;
            default: throw std::invalid_argument("Unsupported element type.");
        }

        // Retrieve global coordinates of nodes
        auto      node_coords_global   = this->node_coords_global(node_coords_system);

        Vec2      best_coords          = initial_guesses[0];
        Precision min_distance_squared = std::numeric_limits<Precision>::max();

        for (const auto& initial_guess : initial_guesses) {
            Precision r = initial_guess(0);
            Precision s = initial_guess(1);

            for (Index iter = 0; iter < max_iter; iter++) {
                auto shape_func         = shape_function(r, s);
                auto shape_deriv        = shape_derivative(r, s);
                auto shape_second_deriv = shape_second_derivative(r, s);

                // Compute position x(r, s) and its derivatives dx/dr, dx/ds, d²x/dr², d²x/ds², d²x/(drds)
                Vec3 x_rs     = Vec3::Zero();
                Vec3 dx_dr    = Vec3::Zero();
                Vec3 dx_ds    = Vec3::Zero();
                Vec3 d2x_dr2  = Vec3::Zero();
                Vec3 d2x_ds2  = Vec3::Zero();
                Vec3 d2x_drds = Vec3::Zero();

                for (Index i = 0; i < N; i++) {
                    x_rs += node_coords_global.row(i) * shape_func(i);
                    dx_dr += node_coords_global.row(i) * shape_deriv(i, 0);
                    dx_ds += node_coords_global.row(i) * shape_deriv(i, 1);
                    d2x_dr2 += node_coords_global.row(i) * shape_second_deriv(i, 0);
                    d2x_ds2 += node_coords_global.row(i) * shape_second_deriv(i, 1);
                    d2x_drds += node_coords_global.row(i) * shape_second_deriv(i, 2);
                }

                // Compute the difference vector
                Vec3 diff = x_rs - p;

                // Compute the derivatives of the squared distance
                Precision dD_dr = diff.dot(dx_dr);
                Precision dD_ds = diff.dot(dx_ds);

                // Compute the Hessian matrix for the squared distance
                Precision d2D_dr2  = dx_dr.dot(dx_dr) + diff.dot(d2x_dr2);
                Precision d2D_ds2  = dx_ds.dot(dx_ds) + diff.dot(d2x_ds2);
                Precision d2D_drds = dx_dr.dot(dx_ds) + diff.dot(d2x_drds);

                // Solve the 2x2 linear system to get Newton-Raphson updates
                StaticMatrix<2, 2> H;
                H << d2D_dr2, d2D_drds, d2D_drds, d2D_ds2;

                Vec2 grad;
                grad << dD_dr, dD_ds;

                Vec2 delta = H.ldlt().solve(-grad);

                // Update r and s
                r += delta(0);
                s += delta(1);

                // Check convergence
                if (delta.norm() < eps) {
                    break;
                }
            }

            if (!clip || in_bounds({r, s})) {
                // Compute the squared distance for this (r, s)
                Vec3      x_rs_final       = this->local_to_global({r, s}, node_coords_system);
                Precision distance_squared = (x_rs_final - p).squaredNorm();

                // Update the best (r, s) if this is the minimum distance so far
                if (distance_squared < min_distance_squared) {
                    min_distance_squared = distance_squared;
                    best_coords          = {r, s};
                }
            }
        }

        if (clip) {
            // check the boundary
            Vec2 best_coords_edge = closest_point_on_boundary(p, node_coords_global);

            // check if its better than the current best
            Precision distance_squared =
                (this->local_to_global(best_coords_edge, node_coords_system) - p).squaredNorm();

            if (distance_squared < min_distance_squared) {
                best_coords = best_coords_edge;
            }
        }

        return best_coords;
    }

    /**
     * @brief Compute the normal vector at a given local coordinate (r, s).
     *
     * @param local Local coordinate in the parametric space of the surface element.
     * @return Vec3 Normal vector at the local coordinate.
     */
    Vec3 normal(const NodeData& node_coords_system, const Vec2& local) const override {
        // Compute the shape function derivatives at the given local coordinates
        auto shape_deriv        = shape_derivative(local(0), local(1));
        auto node_coords_global = this->node_coords_global(node_coords_system);

        // compute the derivative of the current position w.r.t r and s
        Vec3 dx_dr = Vec3::Zero();
        Vec3 dx_ds = Vec3::Zero();

        for (Index i = 0; i < N; i++) {
            dx_dr += node_coords_global.row(i) * shape_deriv(i, 0);
            dx_ds += node_coords_global.row(i) * shape_deriv(i, 1);
        }

        // Compute the normal vector as the cross product of the tangent vectors
        Vec3 normal = dx_dr.cross(dx_ds);
        return normal.normalized();
    }
    /**
     * @brief Computes the integral of each shape function over the global surface area of the element.
     *
     * @param node_coords_system Node coordinates in the global coordinate system.
     * @return StaticMatrix<N, 1> Vector of integrated shape function values over the global surface area.
     *
     * @details The integral of shape functions is computed using numerical integration over the parametric
     *          domain of the surface element, with the Jacobian determinant to map to the global space.
     *          This can be useful for calculating distributed load contributions or surface mass properties.
     */
    DynamicVector shape_function_integral(const NodeData& node_coords_system) const override {
        auto node_coords_global = this->node_coords_global(node_coords_system);

        // Define the integrand for computing the integral of each shape function over the global surface
        std::function<StaticMatrix<N, 1>(Precision, Precision, Precision)> integrand =
            [&](Precision r, Precision s, Precision) -> StaticMatrix<N, 1> {
            // Compute the shape function at the given local coordinates (r, s)
            StaticMatrix<N, 1> shape_func = shape_function(r, s);

            // Compute the 3x2 Jacobian matrix at (r, s)
            auto jac = jacobian(node_coords_global, r, s);

            // Compute the determinant of the Jacobian, which represents the local-to-global area scaling factor
            Precision detJ = (jac.col(0).cross(jac.col(1))).norm();

            // Multiply the shape function values by the determinant to integrate over the global surface area
            return shape_func * detJ;
        };

        // Perform numerical integration using the specified integration scheme
        StaticMatrix<N, 1> shape_function_integral = StaticMatrix<N, 1>::Zero();
        shape_function_integral                    = integration_scheme().integrate(integrand);

        return DynamicVector {shape_function_integral};
    }

    /**
     * @brief finds the closest point on the boundary.
     *
     */
    virtual Vec2 closest_point_on_boundary(const Vec3& global, const StaticMatrix<N, 3>& node_coords) const = 0;

    /**
     * @brief Compute the area of the surface element.
     *
     * @param node_coords_system Node coordinates in the global system.
     * @return Precision Area of the surface element.
     *
     * @details Uses a numerical integration scheme to compute the area of the surface
     *          by integrating the magnitude of the cross product of the tangent vectors
     *          at each integration point.
     */
    Precision area(const NodeData& node_coords_system) const override {
        auto node_coords_global = this->node_coords_global(node_coords_system);

        // Define the integrand for computing surface area
        std::function<Precision(Precision, Precision, Precision)> integrand =
            [&](Precision r, Precision s, Precision) -> Precision {
            // Compute the 3x2 Jacobian
            auto jac = jacobian(node_coords_global, r, s);

            // Compute the cross product of the Jacobian columns
            Vec3 tangent_r = jac.col(0);    // ∂r/∂ξ
            Vec3 tangent_s = jac.col(1);    // ∂r/∂η

            // Compute the magnitude of the cross product: |tangent_r x tangent_s|
            Precision area_element = (tangent_r.cross(tangent_s)).norm();

            return area_element;
        };

        // Perform numerical integration using the specified integration scheme
        Precision total_area = integration_scheme().integrate(integrand);
        return total_area;
    }

    /**
     * @brief Returns a matrix holding the local coordinates of the nodes.
     *
     * @return local coordinates of the nodes.
     */
    virtual StaticMatrix<N, 2> node_coords_local() const = 0;

    /**
     * @brief Extracts the global coordinates using the stored node ids from the global list of node positions
     *
     * @param node_coords_system
     * @return StaticMatrix<N, 3> which holds the global coordinates of the nodes.
     */
    StaticMatrix<N, 3> node_coords_global(const NodeData& node_coords_system) const {
        StaticMatrix<N, 3> res {};
        for (Index i = 0; i < N; i++) {
            for (Index j = 0; j < 3; j++) {
                res(i, j) = node_coords_system(nodeIds[i], j);
            }
        }
        return res;
    }

    /**
     * @brief Get the node IDs associated with the surface element.
     *
     * @return const std::array<ID, N>& Reference to the array of node IDs.
     */
    const std::array<ID, N>& get_node_ids() const {
        return nodeIds;
    }

    /**
     * @brief Provide the integration scheme for the surface element.
     *
     * @return const quadrature::Quadrature& Reference to the quadrature rule used for integration.
     */
    virtual const quadrature::Quadrature& integration_scheme() const = 0;

    /**
     * @brief Get the node IDs associated with the surface element.
     *
     * @return pointer to the node ids
     */
    ID* nodes() override {
        return nodeIds.data();
    }

    /**
     * @brief Apply distributed load to the surface element.
     *
     * @param node_coords Node coordinates in the global system.
     * @param node_loads Node loads in the global system.
     * @param load Vector of the distributed load to be applied.
     */
    void apply_dload(NodeData& node_coords, NodeData& node_loads, Vec3 load) override {
        auto nodal_contributions = shape_function_integral(node_coords);

        // go through the node ids
        for (Index i = 0; i < n_nodes; i++) {
            auto n_id = nodes()[i];
            node_loads(n_id, 0) += nodal_contributions(i) * load(0);
            node_loads(n_id, 1) += nodal_contributions(i) * load(1);
            node_loads(n_id, 2) += nodal_contributions(i) * load(2);
        }
    }

    /**
     * @brief Apply pressure load to the surface element.
     *
     * @param node_coords Node coordinates in the global system.
     * @param node_loads Node loads in the global system.
     * @param load Pressure load to be applied.
     */
    void apply_pload(NodeData& node_coords, NodeData& node_loads, Precision load) override {
        auto node_coords_global = this->node_coords_global(node_coords);

        // Define the integrand for computing the integral of each shape function over the global surface
        std::function<StaticMatrix<N, 3>(Precision, Precision, Precision)> integrand =
            [&](Precision r, Precision s, Precision) -> StaticMatrix<N, 3> {
            // Compute the shape function at the given local coordinates (r, s)
            StaticMatrix<N, 1> shape_func = shape_function(r, s);
            StaticMatrix<N, 3> shape_func_normal{};

            // Compute the 3x2 Jacobian matrix at (r, s)
            auto jac    = jacobian(node_coords_global, r, s);
            Vec3 normal = this->normal(node_coords, {r, s});

            for(Index i = 0; i < N; i++) {
                shape_func_normal(i, 0) = shape_func(i, 0) * normal(0);
                shape_func_normal(i, 1) = shape_func(i, 0) * normal(1);
                shape_func_normal(i, 2) = shape_func(i, 0) * normal(2);
            }

            // Compute the determinant of the Jacobian, which represents the local-to-global area scaling factor
            Precision detJ = (jac.col(0).cross(jac.col(1))).norm();

            // Multiply the shape function values by the determinant to integrate over the global surface area
            return shape_func_normal * load * detJ;
        };

        auto nodal_contributions = integration_scheme().integrate(integrand);

        for(Index local_id = 0; local_id < n_nodes; local_id++) {
            auto n_id = nodes()[local_id];
            node_loads(n_id, 0) += nodal_contributions(local_id, 0);
            node_loads(n_id, 1) += nodal_contributions(local_id, 1);
            node_loads(n_id, 2) += nodal_contributions(local_id, 2);
        }
//        for (ID node_id : *this) {
//            node_loads(node_id, 0) += nodal_contributions(node_coords)(0);
//            node_loads(node_id, 1) += nodal_contributions(node_coords)(1);
//            node_loads(node_id, 2) += nodal_contributions(node_coords)(2);
//        }
    }
};
}    // namespace fem::model
