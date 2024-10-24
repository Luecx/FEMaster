#include "eigen.h"
#include "device.h"
#include "method.h"
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymGEigsShiftSolver.h>
#include <Spectra/MatOp/SymShiftInvert.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

namespace fem::solver {

std::pair<DynamicVector, DynamicMatrix> compute_eigenvalues(SolverDevice device,
                                                            SparseMatrix& mat,
                                                            int num_eigenvalues,
                                                            bool return_vectors) {
    runtime_assert(mat.rows() == mat.cols(), "Matrix must be square");

    const auto N = mat.cols();
    DynamicVector eigenvalues{num_eigenvalues};
    DynamicMatrix eigenvectors{N, num_eigenvalues};
    eigenvalues.setZero();
    eigenvectors.setZero();

#ifndef SUPPORT_GPU
    logging::info(device != CPU, "This build does not support GPU-accelerated eigenvalue computation, falling back to CPU");
    device = CPU;
#else
#endif
    logging::info(device != CPU, "The eigenvalue solver does not support GPU-accelerated eigenvalue computation, falling back to CPU");
    device = CPU;

    logging::info(true, "");
    logging::info(true, "Computing ", num_eigenvalues, " eigenvalues with N=", N);

    logging::up();

#ifdef SUPPORT_GPU
    if (device == GPU) {
        // GPU eigenvalue computation using cuSOLVER (same as before)
        // ...
    }
#endif

    if (device == CPU) {
         // Start the time
         Timer t {};
         t.start();

         logging::info(true, "Using Spectra's SymEigsShiftSolver");

         Spectra::SparseSymShiftSolve<Precision> op(mat);
         Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<Precision>> eigs(op, num_eigenvalues, 2 * num_eigenvalues, 0.0);

         eigs.init();
         eigs.compute(Spectra::SortRule::SmallestMagn);

         if (eigs.info() == Spectra::CompInfo::Successful) {
             // Extract the eigenvalues
             eigenvalues = eigs.eigenvalues();

             // If eigenvectors are requested, extract them
             if (return_vectors) {
                 eigenvectors = eigs.eigenvectors();
             }
         } else {
             logging::error(true, "Spectra eigenvalue computation failed");
         }

         t.stop();
         logging::info(true, "Eigenvalue computation finished on CPU");
         logging::info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
    }

    logging::down();

    // Return the eigenvalues and eigenvectors if requested
    return {eigenvalues, return_vectors ? eigenvectors : DynamicMatrix{}};
}


std::pair<DynamicVector, DynamicMatrix> compute_eigenvalues(SolverDevice device,
                                                                SparseMatrix& A,
                                                                const SparseMatrix& B,
                                                                int num_eigenvalues,
                                                                bool return_vectors) {
    runtime_assert(A.rows() == A.cols() && B.rows() == B.cols() && A.rows() == B.rows(),
                   "Matrices A and B must be square and of the same dimension");

    const auto N = A.cols();
    DynamicVector eigenvalues(num_eigenvalues);
    DynamicMatrix eigenvectors(N, num_eigenvalues);
    eigenvalues.setZero();
    eigenvectors.setZero();

#ifndef SUPPORT_GPU
    logging::info(device != CPU, "This build does not support GPU-accelerated eigenvalue computation, falling back to CPU");
    device = CPU;
#else
#endif
    logging::info(device != CPU, "The eigenvalue solver does not support GPU-accelerated eigenvalue computation, falling back to CPU");
    device = CPU;

    logging::info(true, "");
    logging::info(true, "Computing ", num_eigenvalues, " generalized eigenvalues with N=", N);

    logging::up();

#ifdef SUPPORT_GPU
    if (device == GPU) {
        // GPU eigenvalue computation using cuSOLVER (same as before)
        // ...
    }
#endif
    if (device == CPU) {
        // Start the time
        Timer t {};
        t.start();

        logging::info(true, "Using Spectra's SymGEigsShiftSolver");

        // Define the operation objects
        using OpType  = Spectra::SymShiftInvert<Precision, Eigen::Sparse, Eigen::Sparse>;
        using BOpType = Spectra::SparseSymMatProd<Precision>;

        OpType op(A, B);
        BOpType Bop(B);

        // Construct the solver with the specified shift
        Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::ShiftInvert>
            geigs(op, Bop, num_eigenvalues, std::min((int)(2 * num_eigenvalues), (int)A.rows()), 0);

        geigs.init();
        int nconv = geigs.compute(Spectra::SortRule::LargestMagn);

        if (nconv < num_eigenvalues) {
            logging::warning(true, "Spectra generalized eigenvalue computation failed to converge for all eigenvalues");
            logging::warning(true, "Number of converged eigenvalues: ", nconv);
        }

        if (geigs.info() == Spectra::CompInfo::Successful) {
            // Extract the eigenvalues
            eigenvalues = geigs.eigenvalues();

            // If eigenvectors are requested, extract them
            if (return_vectors) {
                eigenvectors = geigs.eigenvectors();
            }
        } else {
            logging::error(true, "Spectra generalized eigenvalue computation failed");
        }

        t.stop();
        logging::info(true, "Generalized eigenvalue computation finished on CPU");
        logging::info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
    }

    logging::down();

    // Return the eigenvalues and eigenvectors if requested
    return {eigenvalues, return_vectors ? eigenvectors : DynamicMatrix(0,0)};
}


}  // namespace fem::solver
