/**
 * @file device.h
 * @brief Defines the SolverDevice enum used to specify the computational
 *        device (CPU or GPU) for solving linear systems and eigenvalue problems.
 *
 * @details The SolverDevice enum is utilized throughout the solver framework
 *          to indicate whether computations should be performed on the CPU or GPU.
 *          This allows users to control which hardware resources are used during
 *          the computation process.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *         All rights reserved.
 * @date Created on 27.08.2024
 *
 */

#pragma once

namespace fem::solver {

/**
 * @enum SolverDevice
 * @brief Specifies the computational device to use for solving linear systems
 *        or eigenvalue problems.
 *
 * @details This enum allows users to select whether to perform computations on
 *          the CPU or the GPU (using CUDA for the latter). It is used as a parameter
 *          in the solver functions to determine which hardware device will be used
 *          for the computations.
 */
enum SolverDevice {
    GPU, ///< Perform computations on the GPU using CUDA.
    CPU  ///< Perform computations on the CPU.
};

}  // namespace solver
