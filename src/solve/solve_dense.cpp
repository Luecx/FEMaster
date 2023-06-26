////
//// Created by Luecx on 09.06.2023.
////
//#include "solver.h"
//
//namespace solver{
//
//linalg::DenseNMatrix<float>
//    solve(array::Device device, linalg::DenseNMatrix<float>& mat, linalg::DenseNMatrix<float>& rhs) {
//
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
//#ifdef SUPPORT_GPU
//    if (device == array::GPU) {
//        // init gpu if not init yet
//        cuda::manager.create_cuda();
//
//        sol = rhs;
//
//        mat >> array::GPU;
//        sol >> array::GPU;
//
//        // compute workspace size
//        int workspace_size = 0;
//        cusolverDnSgetrf_bufferSize(cuda::manager.handle_cusolver_dn,
//                                    N,
//                                    N,
//                                    mat.address(array::GPU),
//                                    N,
//                                    &workspace_size);
//
//        // additional arrays required
//        array::DynamicSyncedArray<int>   gpu_pivot {N};
//        array::DynamicSyncedArray<float> gpu_works {static_cast<size_t>(workspace_size)};
//        array::DynamicSyncedArray<int>   success {1};
//
//        gpu_pivot.malloc(array::GPU);
//        gpu_works.malloc(array::GPU);
//        success  .malloc(array::BOTH);
//
//        // compute LU factorization
//        runtime_check_cuda(cusolverDnSgetrf(cuda::manager.handle_cusolver_dn,
//                                            N,
//                                            N,
//                                            mat      .address(array::GPU),
//                                            N,
//                                            gpu_works.address(array::GPU),
//                                            gpu_pivot.address(array::GPU),
//                                            success  .address(array::GPU)));
//        success >> array::CPU;
//        log_error(success[0] == 0, "decomposing system not possible");
//
//        // compute solution
//        runtime_check_cuda(cusolverDnSgetrs(cuda::manager.handle_cusolver_dn,
//                                            CUBLAS_OP_N,
//                                            N,
//                                            1,
//                                            mat      .address(array::GPU),
//                                            N,
//                                            gpu_pivot.address(array::GPU),
//                                            sol      .address(array::GPU),
//                                            N,
//                                            success  .address(array::GPU)));
//        sol     >> array::CPU;
//        success >> array::CPU;
//        log_error(success[0] == 0, "solving system not possible");
//
//
//        sol.free(array::GPU);
//        mat.free(array::GPU);
//    }
//
//#endif
//
//    return sol;
//    // TODO
//}
//}
