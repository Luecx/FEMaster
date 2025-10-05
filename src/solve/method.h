/**
* @file method.h
 * @brief Defines the SolverMethod enum used to specify the solver method
 *        (direct or iterative) for solving linear systems and eigenvalue problems.
 *
 * @details The SolverMethod enum allows users to select between direct solvers,
 *          which use methods such as LU decomposition, and iterative solvers,
 *          which use algorithms like Conjugate Gradient or GMRES.
 *
 * @date Created on 27.08.2024
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *         All rights reserved.
 *
 */

#pragma once

namespace fem::solver {

/**
 * @enum SolverMethod
 * @brief Specifies the method to use for solving linear systems or eigenvalue problems.
 *
 * @details This enum allows users to choose between direct methods, which solve the
 *          system using exact methods like LU or QR decomposition, and iterative methods,
 *          which use iterative algorithms like Conjugate Gradient or Lanczos.
 */
enum SolverMethod {
    DIRECT,   ///< Use a direct solver (e.g., LU, QR decomposition).
    INDIRECT  ///< Use an iterative solver (e.g., Conjugate Gradient, GMRES, Lanczos).
};

}  // namespace solver
