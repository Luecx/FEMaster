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
#include <Argparse/argparse.h>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>

int main(int argc, char* argv[]) {
    // Set up argparse
    argparse::ArgumentParser program("FEM Solver");

    program.add_argument("input_file")
        .help("Path to the input file");

    program.add_argument("--ncpus")
        .default_value(1)
        .scan<'i', int>()
        .help("Number of CPUs to use (default: 1)");

    // Parse arguments
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    // Access parsed arguments
    std::string input_file = program.get<std::string>("input_file");
    int ncpus = program.get<int>("--ncpus");

    std::cout << "Input file: " << input_file << std::endl;
    std::cout << "CPU(s)    : " << ncpus << std::endl;

    // store number of cpus in config
    global_config.max_threads = ncpus;


    fem::reader::Reader reader{input_file};
    reader.read();

    return 0;
}
