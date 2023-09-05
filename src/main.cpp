
#include "core/core.h"
#include "cuda/cuda.h"
#include "material/material.h"
#include "model/c3d20.h"
#include "model/c3d8.h"
#include "model/element.h"
#include "model/model.h"
#include "model/quadrature.h"
#include "reader/line.h"
#include "reader/reader.h"
#include "reader/file.h"
#include "solve/solver.h"
#include "loadcase/linear_static.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <chrono>
#include <iomanip>
#include <iostream>

int main(int argc, char* argv[]) {

//    fem::reader::Reader reader{"test1.txt"};
//    reader.read();

    using namespace fem::model;
    using namespace fem::reader;

    const int N = 20;
    const int dofs_per_node = 3;
    const int nodes_per_elem = 8;
    Model model{N*N*N, N*N*N + 4};

    // Generate node coordinates.
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                int n = x * N * N + y * N + z;
                model.set_node(n, x,y,z);
            }
        }
    }

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


    model.activate_material("mat1");

    model.active_material().set_elasticity<fem::material::IsotropicElasticity>(210000, 0);
    model.active_material().set_density(8100);

    model.solid_section("EALL", "mat1");

    model.add_support(0, StaticVector<3>(0,0,0));
    model.add_support(2, StaticVector<3>(0,0,0));
    model.add_support(4, StaticVector<3>(0,0,0));
    model.add_support(6, StaticVector<3>(0,0,0));
    model.add_support(53, StaticVector<3>(0,0,0));
    model.add_cload("NALL", StaticVector<3>(2,0,0));
    model.add_support(30, StaticVector<3>(0,0,0));

    fem::loadcase::LinearStatic lc1{&model};
    lc1.device = solver::CPU;
    lc1.method = solver::DIRECT;
    lc1.run();
    fem::loadcase::LinearStatic lc2{&model};
    lc2.device = solver::CPU;
    lc2.method = solver::INDIRECT;
    lc2.run();
    fem::loadcase::LinearStatic lc3{&model};
    lc3.device = solver::GPU;
    lc3.method = solver::DIRECT;
    lc3.run();
    fem::loadcase::LinearStatic lc4{&model};
    lc4.device = solver::GPU;
    lc4.method = solver::INDIRECT;
    lc4.run();
//
//    NodeData stress;
//    NodeData strain;

//    std::tie(stress, strain) = model.compute_stress_strain(disp);
    return 0;
}
