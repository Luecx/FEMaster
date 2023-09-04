
#include "core/core.h"
#include "cuda/cuda.h"
#include "material/material.h"
#include "model/c3d20.h"
#include "model/c3d8.h"
#include "model/element.h"
#include "model/model.h"
#include "model/quadrature.h"
#include "reader/line.h"
#include "reader/file.h"
#include "solve/solver.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <chrono>
#include <iomanip>
#include <iostream>

int main(int argc, char* argv[]) {


    using namespace fem::model;
    using namespace fem::reader;
//
//    File file{"test1.txt"};
//    while(!file.is_eof()){
//        std::cout << file.next().line() << std::endl;
//    }

//    Line line;
//    line = "* command, key_1 = 3";
//    std:ine line;
//    line = "* command, key_1 = 3";
//    std::cout << line << std::endl;

//    Model model{10,10};
//    model.set_node(0,1,2,3);
//    model.set_element<C3D8>(0,(ID)1,(ID)2,(ID)3,(ID)4,(ID)5,(ID)6,(ID)7);

    const int N = 2;
    const int dofs_per_node = 3;
    const int nodes_per_elem = 8;
    Model model{N*N*N, N*N*N + 4};

    // Generate node coordinates.
    std::cout << "begin setting node coords" << std::endl;
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                int n = x * N * N + y * N + z;
                model.set_node(n, x,y,z);
            }
        }
    }
    std::cout << "finished setting node coords" << std::endl;

    // Material setup.

    int K = N * N * dofs_per_node;
    int M = N * N * N * dofs_per_node - K;
    SparseMatrixBuilder tripplets{};

    // Generate elements and compute stiffness.
    int E = 0;

    SparseMatrix matrix{M, M};

    for (int x = 0; x < N-1; x++) {
        for (int y = 0; y < N-1; y++) {
            for (int z = 0; z < N-1; z++) {
                int e = x * (N-1) * (N-1) + y * (N-1) + z;

                model.set_element<C3D8>(e,
                        (ID)(x)*N*N+(y)*N+(z),
                        (ID)(x+1)*N*N+(y)*N+(z),
                        (ID)(x+1)*N*N+(y+1)*N+(z),
                        (ID)(x)*N*N+(y+1)*N+(z),
                        (ID)(x)*N*N+(y)*N+(z+1),
                        (ID)(x+1)*N*N+(y)*N+(z+1),
                        (ID)(x+1)*N*N+(y+1)*N+(z+1),
                        (ID)(x)*N*N+(y+1)*N+(z+1)
                );
            }
        }
    }

    std::cout << "generated elements" << std::endl;

    model.activate_material("mat1");

    model.active_material().set_elasticity<fem::material::IsotropicElasticity>(210000, 0);
    model.active_material().set_density(8100);

    model.solid_section("EALL", "mat1");

    std::cout << model.node_coords << std::endl;

    std::cout << "added material" << std::endl;
//    model.add_support(0, StaticVector<3>(0,0,0));
    model.add_support(0, StaticVector<3>(0,0,0));
    model.add_support(2, StaticVector<3>(0,0,0));
    model.add_support(4, StaticVector<3>(0,0,0));
    model.add_support(6, StaticVector<3>(0,0,0));
//    model.add_support(53, StaticVector<3>(0,0,0));
    model.add_cload("NALL", StaticVector<3>(2,0,0));
//    model.add_support(30, StaticVector<3>(0,0,0));

    std::cout << "added material and loads" << std::endl;


    Timer timer;  // Create an instance of Timer

    timer.start();
    auto unconstrained = model.build_unconstrained_index_matrix();
    timer.stop();
    std::cout << "Time taken for build_unconstrained_index_matrix: " << timer.elapsed() << " ms" << std::endl;

    timer.start();
    auto supp_vec = model.build_support_vector(unconstrained);
    timer.stop();
    std::cout << "Time taken for build_support_vector: " << timer.elapsed() << " ms" << std::endl;

    timer.start();
    auto load_vec = model.build_load_vector(unconstrained);
    timer.stop();
    std::cout << "Time taken for build_load_vector: " << timer.elapsed() << " ms" << std::endl;

    timer.start();
    auto stiffness = model.build_stiffness_matrix(unconstrained);
    timer.stop();
    std::cout << "Time taken for build_stiffness_matrix: " << timer.elapsed() << " ms" << std::endl;

    timer.start();
    auto impl_load_vec = model.build_implicit_load_vector(stiffness, supp_vec);
    timer.stop();
    std::cout << "Time taken for build_implicit_load_vector: " << timer.elapsed() << " ms" << std::endl;

    timer.start();
    auto constrained = model.build_constrained_index_matrix(supp_vec);
    timer.stop();
    std::cout << "Time taken for build_constrained_index_matrix: " << timer.elapsed() << " ms" << std::endl;

    timer.start();
    auto mapping_vec = model.build_mapping_vector(unconstrained, constrained);
    timer.stop();
    std::cout << "Time taken for build_mapping_vector: " << timer.elapsed() << " ms" << std::endl;

    load_vec += impl_load_vec;

    timer.start();
    auto reduced_stiffness = model.build_reduced_stiffness(mapping_vec, stiffness);
    timer.stop();
    std::cout << "Time taken for build_reduced_stiffness: " << timer.elapsed() << " ms" << std::endl;

    timer.start();
    auto reduced_load = model.build_reduced_load(mapping_vec, load_vec);
    timer.stop();
    std::cout << "Time taken for build_reduced_load: " << timer.elapsed() << " ms" << std::endl;

    auto res1 = solver::solve_iter(solver::CPU, reduced_stiffness, reduced_load);
    auto disp = model.build_global_displacement(constrained, res1);

    std::cout << disp << std::endl;

    NodeData stress;
    NodeData strain;

    std::tie(stress, strain) = model.compute_stress_strain(disp);
    std::cout << stress << std::endl;
    std::cout << strain << std::endl;
    ElementData compliance = model.compute_compliance(disp);
    std::cout << compliance << std::endl;

    return 0;
}
