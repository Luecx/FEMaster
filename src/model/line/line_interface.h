/******************************************************************************
* @file line_interface.h
 * @brief Interface definition for isoparametric line elements in finite element models.
 *
 * @details This file provides a generic interface for line elements with various node configurations,
 *          such as linear, quadratic, and cubic line elements (2, 3, or 4 nodes). The interface defines
 *          essential functions for computing shape functions, their derivatives, and coordinate transformations
 *          between local and global coordinates. Additionally, utility functions such as length computation
 *          and global-to-local coordinate mapping are provided to facilitate geometric analysis of line elements.
 *
 * @tparam N  Number of nodes in the line element (must be 2, 3, or 4).
 * @tparam CR Coordinate range type: MINUS_ONE_TO_ONE (default) or ZERO_TO_ONE.
 *
 * @note This file is part of the finite element model library and is used to build higher-level element
 *       types and analysis routines. Users can create custom line elements by extending this interface.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *         All rights reserved.
 * @date Created on 30.09.2024
 *
 ******************************************************************************/

#pragma once

#include <array>
#include <Eigen/Dense>
#include <functional>

#include "../element.h"
#include "../element_solid.h"
#include "../../math/quadrature.h"
#include "../../core/types.h"

namespace fem::model {

/**
 * @enum IsoParametricLineRange
 * @brief Specifies the coordinate range for isoparametric line elements.
 *
 * @details This enum is used to select the coordinate range type for isoparametric line elements,
 *          either from -1 to 1 or from 0 to 1. This affects the shape functions and integration
 *          domains used in computations involving line elements.
 */
enum IsoParametricLineRange {
    MINUS_ONE_TO_ONE, ///< Coordinate range from -1 to 1.
    ZERO_TO_ONE       ///< Coordinate range from 0 to 1.
};

/**
 * @brief Interface for isoparametric line elements in finite element analysis.
 *
 * @tparam N  Number of nodes in the line element (must be 2, 3, or 4).
 * @tparam CR Coordinate range type: MINUS_ONE_TO_ONE (default) or ZERO_TO_ONE.
 *
 * @details This interface defines the essential functions and data for line elements used in finite element models.
 *          It provides methods for computing shape functions, derivatives, element length, and coordinate transformations
 *          between local (parametric) and global coordinate systems.
 */
template<Index N, IsoParametricLineRange CR = MINUS_ONE_TO_ONE>
struct LineInterface {
    /**
     * @brief Number of nodes in the line element.
     */
    constexpr static Index num_nodes = N;

    static_assert(N == 2 || N == 3 || N == 4, "Number of nodes (N) must be 2, 3, or 4.");

    /**
     * @brief IDs of the nodes that make up the line element.
     */
    std::array<ID, N> node_ids;

    /**
     * @brief Constructor for LineInterface.
     *
     * @param pNodeIds Array of node IDs for the line element.
     */
    LineInterface(const std::array<ID, N>& pNodeIds) : node_ids(pNodeIds) {}

    /**
     * @brief Compute the shape functions at a given local coordinate.
     *
     * @param r Local coordinate (parametric coordinate along the line element).
     * @return Eigen::Matrix<Precision, N, 1> Vector of shape function values at coordinate r.
     */
    virtual Eigen::Matrix<Precision, N, 1> shape_function(Precision r) const = 0;

    /**
     * @brief Compute the first derivative of shape functions at a given local coordinate.
     *
     * @param r Local coordinate (parametric coordinate along the line element).
     * @return Eigen::Matrix<Precision, N, 1> Vector of first derivatives of shape functions at coordinate r.
     */
    virtual Eigen::Matrix<Precision, N, 1> shape_derivative(Precision r) const = 0;

    /**
     * @brief Compute the second derivative of shape functions at a given local coordinate.
     *
     * @param r Local coordinate (parametric coordinate along the line element).
     * @return Eigen::Matrix<Precision, N, 1> Vector of second derivatives of shape functions at coordinate r.
     */
    virtual Eigen::Matrix<Precision, N, 1> shape_second_derivative(Precision r) const = 0;

    /**
     * @brief Retrieve the global coordinates of the nodes associated with the line element.
     *
     * @param node_coords_system System-wide node coordinates data structure.
     * @return StaticMatrix<N, 3> Matrix containing the global coordinates of the element's nodes.
     *
     * @details This function maps the node IDs of the line element to their corresponding global coordinates
     *          in the provided node coordinate system. It returns an N x 3 matrix where each row corresponds
     *          to a node and columns represent the x, y, and z coordinates.
     */
    virtual StaticMatrix<N, 3> node_coords_global(const NodeData& node_coords_system) const {
        StaticMatrix<N, 3> res {};
        for (Index i = 0; i < N; i++) {
            for(Index j = 0; j < 3; j++) {
                res(i, j) = node_coords_system(node_ids[i], j);
            }
        }
        return res;
    }

    /**
     * @brief Compute the length of the line element using numerical integration.
     *
     * @param node_coords_system Node coordinates in the global coordinate system.
     * @return Precision Length of the line element.
     *
     * @details This function calculates the length of the line element by integrating over the element's
     *          parametric domain using a quadrature rule. The length is computed by integrating the norm
     *          of the derivative of the position vector with respect to the local coordinate r.
     */
    Precision length(const NodeData& node_coords_system) const {
        using namespace fem::quadrature;

        // Select the appropriate domain
        Domain domain = (CR == MINUS_ONE_TO_ONE) ? DOMAIN_ISO_LINE_A : DOMAIN_ISO_LINE_B;

        // Create the quadrature rule (order 2, 3-point quadrature for length computation)
        Quadrature quadrature(domain, ORDER_QUADRATIC);

        // Get node global coordinates
        StaticMatrix<N, 3> node_coords_global = this->node_coords_global(node_coords_system);

        // Integrate to compute the length
        return quadrature.integrate([&](Precision r, Precision, Precision) -> Precision {
            StaticVector<3> dx_dr = StaticVector<3>::Zero();
            Eigen::Matrix<Precision, N, 1> dN_dr = this->shape_derivative(r);

            for (Index i = 0; i < N; i++) {
                dx_dr += dN_dr(i) * node_coords_global.row(i);
            }
            return dx_dr.norm();
        });
    }

    /**
     * @brief Find the local coordinate 'r' corresponding to the closest point on the line element to a given global point.
     *
     * @param p The point in global coordinates to which the closest point on the line element is sought.
     * @param node_coords_system Node coordinates in the global coordinate system.
     * @return Precision The local coordinate 'r' corresponding to the closest point on the line element to point 'p'.
     *
     * @details This function uses an iterative method (Newton-Raphson) to find the parametric coordinate 'r' that
     *          minimizes the distance between the point 'p' and the position vector of the line element. Multiple
     *          initial guesses are used to ensure convergence to the global minimum distance.
     */
    Precision global_to_local(const StaticVector<3>& p, const NodeData& node_coords_system, bool clip=true) const {
        constexpr int max_iter = 32;
        constexpr Precision eps = 1e-12;
        std::vector<Precision> initial_guesses;

        // Determine initial guesses based on N
        if (N == 3) {
            initial_guesses = { min_r(), max_r() };
        } else if (N > 3) {
            initial_guesses = { min_r(), (min_r() + max_r()) / 2.0, max_r() };
        } else {
            initial_guesses = { (min_r() + max_r()) / 2.0 };
        }

        // Retrieve global coordinates of nodes
        auto node_coords_global = this->node_coords_global(node_coords_system);

        Precision best_r = min_r();
        Precision min_distance_squared = std::numeric_limits<Precision>::max();

        for (auto initial_r : initial_guesses) {

            Precision r = initial_r;

            for (Index iter = 0; iter < max_iter; iter++) {
                // Compute shape functions and derivatives at current r
                Eigen::Matrix<Precision, N, 1> N_vals = shape_function(r);
                Eigen::Matrix<Precision, N, 1> dN_dr = shape_derivative(r);
                Eigen::Matrix<Precision, N, 1> d2N_dr2 = shape_second_derivative(r);

                // Compute position x(r) and its derivatives dx/dr and d²x/dr²
                StaticVector<3> x_r = StaticVector<3>::Zero();
                StaticVector<3> dx_dr = StaticVector<3>::Zero();
                StaticVector<3> d2x_dr2 = StaticVector<3>::Zero();

                for (Index i = 0; i < N; i++) {
                    x_r += N_vals(i) * node_coords_global.row(i);
                    dx_dr += dN_dr(i) * node_coords_global.row(i);
                    d2x_dr2 += d2N_dr2(i) * node_coords_global.row(i);
                }

                // Compute the difference vector
                StaticVector<3> diff = x_r - p;

                // Compute the derivative of the squared distance
                Precision dD_dr = diff.dot(dx_dr);

                // Compute the second derivative of the squared distance
                Precision d2D_dr2 = (dx_dr.dot(dx_dr) + diff.dot(d2x_dr2));

                // Newton-Raphson update
                Precision delta_r = -dD_dr / d2D_dr2;

                // Update r
                r += delta_r;

                if (clip) {
                    // Ensure r stays within bounds
                    if (r < min_r()) r = min_r();
                    if (r > max_r()) r = max_r();
                }

                // Check convergence
                if (std::abs(delta_r) < eps) {
                    break;
                }
            }

            // Compute the squared distance for this r
            StaticVector<3> x_r = this->local_to_global(r, node_coords_system);
            Precision distance_squared = (x_r - p).squaredNorm();

            // Update the best r if this is the minimum distance so far and within bounds
            if (distance_squared < min_distance_squared) {
                min_distance_squared = distance_squared;
                best_r = r;
            }
        }

        return best_r;
    }

    /**
     * @brief Map a local coordinate 'r' to the corresponding global coordinate on the line element.
     *
     * @param r Local coordinate along the line element.
     * @param node_coords_system Node coordinates in the global coordinate system.
     * @return StaticVector<3> Global coordinate corresponding to the local coordinate 'r'.
     *
     * @details This function computes the global position vector at a given local coordinate 'r' by evaluating
     *          the shape functions at 'r' and summing the contributions from each node.
     */
    StaticVector<3> local_to_global(Precision r, const NodeData& node_coords_system) const {
        auto node_coords_global = this->node_coords_global(node_coords_system);
        StaticVector<3> res = StaticVector<3>::Zero();
        Eigen::Matrix<Precision, N, 1> N_vals = shape_function(r);
        for (Index i = 0; i < N; i++) {
            res += N_vals(i) * node_coords_global.row(i);
        }
        return res;
    }

    /**
     * @brief Get the minimum value of the local coordinate 'r' for the line element.
     *
     * @return Precision Minimum value of 'r' based on the coordinate range type.
     */
    constexpr Precision min_r() const { return CR == MINUS_ONE_TO_ONE ? -1 : 0; }

    /**
     * @brief Get the maximum value of the local coordinate 'r' for the line element.
     *
     * @return Precision Maximum value of 'r' (always 1).
     */
    constexpr Precision max_r() const { return CR == MINUS_ONE_TO_ONE ?  1 : 1; }

    /**
     * @brief Get the node IDs associated with the line element.
     *
     * @return const std::array<ID, N>& Reference to the array of node IDs.
     */
    const std::array<ID, N>& get_node_ids() const {
        return node_ids;
    }
};
}