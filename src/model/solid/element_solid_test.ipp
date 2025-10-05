/**
 * @file element_solid.ipp
 * @brief Implementation of the SolidElement class template. This file contains
 * the definitions for methods declared in the SolidElement class, including
 * the computation of the strain-displacement matrix, Jacobian, stiffness and
 * mass matrices, and other element-level calculations.
 *
 * @date Created on 12.06.2023
 * @author Finn Eggers
 */

#pragma once

namespace fem::model {

template<Index N>
template<class ElementType>
bool
SolidElement<N>::test_implementation(bool print) {
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
    if (print)
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
                    if (print)
                        std::cout << "Sum of shape functions at (r, s, t) = (" << r << ", " << s << ", " << t
                                  << ") is not 1. Actual sum: " << sum << std::endl;
                    return false;
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
                            if (print)
                                std::cout << "Mismatch in derivative at (r, s, t) = (" << r << ", " << s << ", " << t
                                          << ") in direction " << d
                                          << " from shape function " << j
                                          << ". True derivative: " << true_derivatives(j, d)
                                          << ", Finite Difference: " << finite_diff_derivatives(j, d) << std::endl;
                            return false;
                        }
                    }
                }
            }
        }
    }
    return true;
}

}  // namespace fem::model
