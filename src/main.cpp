
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

    fem::reader::Reader reader{std::string(argv[1])};
    reader.read();

    return 0;
}
