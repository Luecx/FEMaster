
#include "core/core.h"
#include "cuda/cuda.h"
#include "loadcase/linear_static.h"
#include "material/material.h"
#include "math/interpolate.h"
#include "math/quadrature.h"
#include "model/c3d10.h"
#include "model/c3d15.h"
#include "model/c3d20.h"
#include "model/c3d4.h"
#include "model/c3d6.h"
#include "model/c3d8.h"
#include "model/element.h"
#include "model/model.h"
#include "reader/file.h"
#include "reader/line.h"
#include "reader/reader.h"
#include "solve/solver.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>

int main(int argc, char* argv[]) {




//    const int N = 16;

//    NodeData xyz(N, 3);
//    NodeData values(N, 1);
//
//    for (int i = 0; i < N; ++i) {
//        xyz(i, 0) = 0.5 + 1 / (2 * sqrt(3)) - 1 / (sqrt(3)) * ((i / 4) % 2);
//        xyz(i, 1) = 0.5 + 1 / (2 * sqrt(3)) - 1 / (sqrt(3)) * ((i / 2) % 2);
//        xyz(i, 2) = 1.5 + 1 / (2 * sqrt(3)) - 1 / (sqrt(3)) * (i % 2) - (i / 8) % 2;
//        values(i, 0) = 400;
//    }
//    std::cout << xyz << std::endl;
//    std::cout << values << std::endl;
//    DynamicVector r2 {1};
//    DynamicVector res = fem::math::interpolate::interpolate<fem::math::interpolate::CONSTANT>(xyz, values, {0,0,0}, &r2);
//    std::cout << res << std::endl;
//    std::cout << r2 << std::endl;

    fem::reader::Reader reader{std::string(argv[1])};
    reader.read();

    return 0;
}
