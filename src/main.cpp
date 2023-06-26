
#include "core/core.h"
#include "cuda/cuda.h"
#include "material/material.h"
#include "model/c3d20.h"
#include "model/c3d8.h"
#include "model/element.h"
#include "model/model.h"
#include "model/quadrature.h"
#include "solve/solver.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <chrono>
#include <iomanip>
#include <iostream>

int main(int argc, char* argv[]) {


    using namespace fem::model;

    Model model{10,10};
    model.set_node(0,1,2,3);
    model.set_element<C3D8>(0,(ID)1,(ID)2,(ID)3,(ID)4,(ID)5,(ID)6,(ID)7);

//    const int N = 30;
//    const int dofs_per_node = 3;
//    const int nodes_per_elem = 8;
//    Model model{N*N*N, N*N*N};
//
//    // Generate node coordinates.
//    std::cout << "begin setting node coords" << std::endl;
//    for (int x = 0; x < N; x++) {
//        for (int y = 0; y < N; y++) {
//            for (int z = 0; z < N; z++) {
//                int n = x * N * N + y * N + z;
//                model.node_coords(n,0) = x;
//                model.node_coords(n,1) = y;
//                model.node_coords(n,2) = z;
//            }
//        }
//    }
//    std::cout << "finished setting node coords" << std::endl;
//
//    // Material setup.
//    fem::material::MaterialPtr mat = fem::material::MaterialPtr(new fem::material::Material("my_mat"));
//    mat->set_elasticity<fem::material::IsotropicElasticity>(210000, 0.3);
//    mat->set_density(8100);
//
//
//    int K = N * N * dofs_per_node;
//    int M = N * N * N * dofs_per_node - K;
//    SparseMatrixBuilder tripplets{};
//
//    // Generate elements and compute stiffness.
//    int E = 0;
//
//    SparseMatrix matrix{M, M};
//
//    for (int x = 0; x < N-1; x++) {
//        for (int y = 0; y < N-1; y++) {
//            for (int z = 0; z < N-1; z++) {
//                int e = x * (N-1) * (N-1) + y * (N-1) + z;
//                C3D8 elem(e, &model, {
//                                         (ID)(x)*N*N+(y)*N+(z),
//                                         (ID)(x+1)*N*N+(y)*N+(z),
//                                         (ID)(x+1)*N*N+(y+1)*N+(z),
//                                         (ID)(x)*N*N+(y+1)*N+(z),
//                                         (ID)(x)*N*N+(y)*N+(z+1),
//                                         (ID)(x+1)*N*N+(y)*N+(z+1),
//                                         (ID)(x+1)*N*N+(y+1)*N+(z+1),
//                                         (ID)(x)*N*N+(y+1)*N+(z+1)
//                                     });
//                elem.set_material(mat);
//
//                {
//                    E ++;
//                    if ( E % 1000 == 0){
//                        std::cout << E << std::endl;
//                    }
//                };
//
//                // once we collect 10M tripplets, we merge them
//                if(tripplets.size() > 1e9){
//                    std::cout << "collapsing elements "<< std::endl;
//                    matrix.insertFromTriplets(tripplets.begin(), tripplets.end());
//                    tripplets.clear();
//                }
//
//                // Get stiffness matrix for current element.
//                DynamicMatrix K = elem.stiffness();
//
//                // Assemble into global matrix.
//                for (int i = 0; i < nodes_per_elem; i++) {
//                    for (int j = 0; j < nodes_per_elem; j++) {
//                        for (int idof = 0; idof < dofs_per_node; idof++) {
//                            for (int jdof = 0; jdof < dofs_per_node; jdof++) {
//                                size_t global_i = elem.node_ids[i] * dofs_per_node + idof;
//                                size_t global_j = elem.node_ids[j] * dofs_per_node + jdof;
//
//                                // Assuming global matrix is large enough to accommodate all entries.
//                                if (global_i < M && global_j < M){
//                                    {
//                                        tripplets.emplace_back(global_i, global_j, K(i * dofs_per_node + idof, j * dofs_per_node + jdof));
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    matrix.insertFromTriplets(tripplets.begin(), tripplets.end());
//
//    DynamicVector rhs{matrix.rows()};
//    for(int i = 0; i < rhs.size(); i++){
//        rhs[i] = rand() / (Precision)(RAND_MAX) / 1000.f;
//    }
//
//    matrix.finalize();
//    solver::solve_iter(solver::GPU, matrix, rhs);
//    solver::solve_iter(solver::GPU, matrix, rhs);


    return 0;
}
