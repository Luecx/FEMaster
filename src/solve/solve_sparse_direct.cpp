//
//#include "solver.h"
//
//#ifdef SUPPORT_GPU
//
//namespace solver{
//
//linalg::DenseNMatrix<float>
//    solve(array::Device device, linalg::ImmutableSparseMatrix<float>& mat, linalg::DenseNMatrix<float>& rhs) {
//
//    runtime_assert(mat.m() == mat.n(), "matrix must be square");
//    runtime_assert(rhs.m() == mat.m(), "missmatch of rhs and matrix");
//    runtime_assert(rhs.n() == 1, "can only solve one equation at a time");
//
//    const auto N = rhs.m();
//    linalg::DenseNMatrix<float> sol{N};
//
//#ifndef SUPPORT_GPU
//    log_info(device != array::CPU, "This build does not support gpu-accelerated solving, falling back to cpu");
//    device = array::CPU;
//#endif
//
//    log_info(true, "");
//    log_info(true, "Solving system with N=", N, " nnz=", mat.get_non_zero_values().size(), " using Cholesky decomposition");
//
//#ifdef SUPPORT_GPU
//    if (device == array::GPU) {
//        // init gpu if not init yet
//        cuda::manager.create_cuda();
//
//        // Move all parts of the matrix and vector to the GPU
//        mat.get_non_zero_values() >> array::GPU;
//        mat.get_row_extents()     >> array::GPU;
//        mat.get_col_indices()     >> array::GPU;
//        sol = rhs;
//        sol >> array::GPU;
//
//        // get number of non-zero values
//        int nnz = static_cast<int>(mat.get_non_zero_values().size());
//
//        // create matrix descriptor
//        cusparseMatDescr_t descr;
//        runtime_check_cuda(cusparseCreateMatDescr(&descr));
//        runtime_check_cuda(cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL));
//        runtime_check_cuda(cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO));
//
//        // stop the time
//        Timer t {};
//        t.start();
//        int singularity;
//        runtime_check_cuda(cusolverSpScsrlsvchol(cuda::manager.handle_cusolver_sp,
//                                                 N,
//                                                 nnz,
//                                                 descr,
//                                                 mat.get_non_zero_values().address(array::GPU),
//                                                 mat.get_row_extents().address(array::GPU),
//                                                 mat.get_col_indices().address(array::GPU),
//                                                 sol.address(array::GPU),
//                                                 0,    // tolerance
//                                                 3,    // reorder
//                                                 sol.address(array::GPU),
//                                                 &singularity));
//        t.stop();
//        log_error(singularity == -1, "decomposing system not possible");
//        log_info(true, "Running PCG method finished");
//        log_info(true, "Elapsed time: ", t.elapsed()," ms");
//
//        sol >> array::CPU;
//
//        rhs                      .free(array::GPU);
//        sol                      .free(array::GPU);
//        mat.get_non_zero_values().free(array::GPU);
//        mat.get_row_extents()    .free(array::GPU);
//        mat.get_col_indices()    .free(array::GPU);
//
//        // destroy matrix descriptor
//        runtime_check_cuda(cusparseDestroyMatDescr(descr));
//    }
//#endif
//
//    //...
//    if (device == array::CPU) {
//        // start the time
//        Timer t {};
//        t.start();
//
//        // First, we need to convert your ImmutableSparseMatrix into an Eigen SparseMatrix
//        Eigen::SparseMatrix<float> eigen_mat(N, N);
//
//        // Create a triplet list for easy construction of the sparse matrix
//        std::vector<Eigen::Triplet<float>> tripletList;
//        tripletList.reserve(mat.get_non_zero_values().size());
//
//        for (size_t i = 0; i < N; ++i) {
//            auto row_start = mat.get_row_extents()[i];
//            auto row_end   = mat.get_row_extents()[i + 1];
//            for (auto j = row_start; j < row_end; ++j) {
//                tripletList.emplace_back(i, mat.get_col_indices()[j], mat.get_non_zero_values()[j]);
//            }
//        }
//        eigen_mat.setFromTriplets(tripletList.begin(), tripletList.end());
//
//        // Convert your DenseVector into an Eigen Vector
//        Eigen::VectorXf eigen_rhs(N);
//        for (size_t i = 0; i < rhs.size(); ++i) {
//            eigen_rhs(i) = rhs[i];
//        }
//
//        // Now we use Eigen's SimplicialLDLT solver to solve the system
//        Eigen::SimplicialLLT<Eigen::SparseMatrix<float>> solver {eigen_mat};
//        solver.compute(eigen_mat);
//        log_error(solver.info() == Eigen::Success, "Decomposition failed");
//        Eigen::VectorXf eigen_sol = solver.solve(eigen_rhs);
//        log_error(solver.info() == Eigen::Success, "Solving failed");
//
//        t.stop();
//        log_info(true, "Solving finished");
//        log_info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
//
//        // Finally, copy the result back into your rhs vector.
//        for (size_t i = 0; i < rhs.size(); ++i) {
//            sol[i] = eigen_sol(i);
//        }
//    }
//    return sol;
//    //...
//}
//
//}
//
//#endif