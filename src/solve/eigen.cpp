#include "eigen.h"
#include "device.h"
#include "method.h"
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymGEigsShiftSolver.h>
#include <Spectra/MatOp/SymShiftInvert.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Eigen/Eigenvalues>

namespace fem::solver {

std::pair<DynamicVector, DynamicMatrix> compute_eigenvalues(SolverDevice device,
                                                            SparseMatrix& mat,
                                                            int num_eigenvalues,
                                                            bool return_vectors) {
    runtime_assert(mat.rows() == mat.cols(), "Matrix must be square");

    const auto N = mat.cols();

    if (num_eigenvalues > N) {
        logging::warning(true, "Number of requested eigenvalues is greater than the matrix size, reducing to N=", N);
        num_eigenvalues = N;
    }

    DynamicVector eigenvalues(num_eigenvalues);
    DynamicMatrix eigenvectors(N, num_eigenvalues);
    eigenvalues.setZero();
    eigenvectors.setZero();

#ifndef SUPPORT_GPU
    logging::info(device != CPU, "This build does not support GPU-accelerated eigenvalue computation, falling back to CPU");
    device = CPU;
#endif

    logging::info(true, "Computing ", num_eigenvalues, " eigenvalues with N=", N);
    logging::up();

#ifdef SUPPORT_GPU
    if (device == GPU) {
        // GPU eigenvalue computation using cuSOLVER (same as before)
        // ...
    }
#endif

    if (device == CPU) {
        Timer t{};
        t.start();

        if (num_eigenvalues == N) {
            logging::info(true, "Using Eigen's SelfAdjointEigenSolver for full eigenvalue computation");
            Eigen::SelfAdjointEigenSolver<DynamicMatrix> solver(mat);

            if (solver.info() == Eigen::Success) {
                eigenvalues = solver.eigenvalues();
                if (return_vectors) {
                    eigenvectors = solver.eigenvectors();
                }
            } else {
                logging::error(true, "Eigen full eigenvalue computation failed");
            }
        } else {
            logging::info(true, "Using Spectra's SymEigsShiftSolver for partial eigenvalue computation");

            int vecs = std::min((int)(2 * num_eigenvalues), (int)N);

            Spectra::SparseSymShiftSolve<Precision> op(mat);
            Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<Precision>> eigs(op, num_eigenvalues, vecs, 0.0);

            eigs.init();
            eigs.compute(Spectra::SortRule::SmallestMagn);

            if (eigs.info() == Spectra::CompInfo::Successful) {
                eigenvalues = eigs.eigenvalues();
                if (return_vectors) {
                    eigenvectors = eigs.eigenvectors();
                }
            } else {
                logging::error(true, "Spectra eigenvalue computation failed");
            }
        }

        t.stop();
        logging::info(true, "Eigenvalue computation finished on CPU");
        logging::info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
    }

    logging::down();
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

    if (num_eigenvalues > N) {
        logging::warning(true, "Number of requested eigenvalues is greater than the matrix size, reducing to N=", N);
        num_eigenvalues = N;
    }

    DynamicVector eigenvalues(num_eigenvalues);
    DynamicMatrix eigenvectors(N, num_eigenvalues);
    eigenvalues.setZero();
    eigenvectors.setZero();

#ifndef SUPPORT_GPU
    logging::info(device != CPU, "This build does not support GPU-accelerated eigenvalue computation, falling back to CPU");
    device = CPU;
#endif

    logging::info(true, "Computing ", num_eigenvalues, " generalized eigenvalues with N=", N);
    logging::up();

#ifdef SUPPORT_GPU
    if (device == GPU) {
        // GPU eigenvalue computation using cuSOLVER (same as before)
        // ...
    }
#endif

    if (device == CPU) {
        Timer t{};
        t.start();

        if (num_eigenvalues == N) {
            logging::info(true, "Using Eigen's Generalized SelfAdjointEigenSolver for full generalized eigenvalue computation");
            Eigen::GeneralizedSelfAdjointEigenSolver<DynamicMatrix> solver(A, B);

            if (solver.info() == Eigen::Success) {
                eigenvalues = solver.eigenvalues();
                if (return_vectors) {
                    eigenvectors = solver.eigenvectors();
                }
            } else {
                logging::error(true, "Eigen full generalized eigenvalue computation failed");
            }
        } else {
            logging::info(true, "Using Spectra's SymGEigsShiftSolver for partial generalized eigenvalue computation");

            using OpType = Spectra::SymShiftInvert<Precision, Eigen::Sparse, Eigen::Sparse>;
            using BOpType = Spectra::SparseSymMatProd<Precision>;

            OpType op(A, B);
            BOpType Bop(B);

            Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::ShiftInvert>
                geigs(op, Bop, num_eigenvalues, std::min((int)(2 * num_eigenvalues), (int)A.rows()), 0);

            geigs.init();
            int nconv = geigs.compute(Spectra::SortRule::LargestMagn);

            if (nconv < num_eigenvalues) {
                logging::warning(true, "Spectra generalized eigenvalue computation failed to converge for all eigenvalues");
                logging::warning(true, "Number of converged eigenvalues: ", nconv);
            }

            if (geigs.info() == Spectra::CompInfo::Successful) {
                eigenvalues = geigs.eigenvalues();
                if (return_vectors) {
                    eigenvectors = geigs.eigenvectors();
                }
            } else {
                logging::error(true, "Spectra generalized eigenvalue computation failed");
            }
        }

        t.stop();
        logging::info(true, "Generalized eigenvalue computation finished on CPU");
        logging::info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
    }

    logging::down();
    return {eigenvalues, return_vectors ? eigenvectors : DynamicMatrix{}};
}

}  // namespace fem::solver
