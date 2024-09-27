// #include "core/core.h"
// #include "core/logging.h"
// #include "cuda/cuda.h"
// #include "loadcase/linear_static.h"
// #include "material/material.h"
// #include "math/interpolate.h"
// #include "math/quadrature.h"
// #include "model/c3d10.h"
// #include "model/c3d15.h"
// #include "model/c3d20.h"
// #include "model/c3d4.h"
// #include "model/c3d6.h"
// #include "model/c3d8.h"
// #include "model/element.h"
// #include "model/model.h"
// #include "reader/file.h"
// #include "reader/line.h"
// #include "reader/reader.h"
// #include "solve/solver.h"

#include "model/surface/surface.h"
#include "model/surface/surface3.h"

int main(int argc, char* argv[]) {
    NodeData node_coords{3, 3};
    node_coords << 0.0, 0.0, 1.0,
                    4.0, 0.0, 0.0,
                    0.0, 1.0, 0.0;
    fem::model::Surface3 surface{{0, 1, 2}};

    auto point = StaticVector<3>{3,4,2};
    auto local_coords1 = surface.global_to_local({3,4,2}, node_coords, false);
    auto local_coords2 = surface.global_to_local({3,4,2}, node_coords, true);

std::cout << local_coords1 << std::endl;
    std::cout << local_coords2 << std::endl;

    // auto projected = surface.local_to_global(local_coords, node_coords);
    // auto dist = (point - projected).norm();
    //
    // std::cout << dist << std::endl;
    // std::cout << surface.shape_function(local_coords[0], local_coords[1]) << std::endl;
}

// #include <Eigen/Core>
// #include <Eigen/Sparse>
// #include <argparse/argparse.hpp>
// #include <chrono>
// #include <functional>
// #include <iomanip>
// #include <iostream>
// #include <random>
// #include <filesystem>
//
// int main(int argc, char* argv[]) {
//     // Set up argparse
//     argparse::ArgumentParser program("FEM Solver");
//
//     program.add_argument("input_file")
//         .help("Path to the input file (.inp is optional).");
//
//     program.add_argument("--ncpus")
//         .default_value(1)
//         .scan<'i', int>()
//         .help("Number of CPUs to use (default: 1)");
//
//     // Parse arguments
//     try {
//         program.parse_args(argc, argv);
//     } catch (const std::runtime_error& err) {
//         std::cerr << err.what() << std::endl;
//         std::cerr << program;
//         return 1;
//     }
//
//     // Access parsed arguments
//     std::string input_file = program.get<std::string>("input_file");
//     int ncpus = program.get<int>("--ncpus");
//
//     // Ensure the input file has ".inp" extension
//     if (input_file.find(".inp") == std::string::npos) {
//         input_file += ".inp";
//     }
//
//     // Check if input file exists
//     if (!std::filesystem::exists(input_file)) {
//         std::cerr << "Error: Input file '" << input_file << "' does not exist." << std::endl;
//         return 1;
//     }
//
//     // Create the output file by replacing ".inp" with ".res"
//     std::string output_file = input_file;
//     size_t pos = output_file.find(".inp");
//     if (pos != std::string::npos) {
//         output_file.replace(pos, 4, ".res");
//     }
//
//     // Logging input and output file information
//     logging::info(true, "");
//     logging::info(true, "Input file: ", input_file);
//     logging::info(true, "Output file: ", output_file);
//     logging::info(true, "CPU(s)    : ", ncpus);
//     logging::info(true, "");
//
//     // Store number of CPUs in config
//     global_config.max_threads = ncpus;
//
//     // Read the input file using the reader
//     fem::reader::Reader reader{input_file, output_file};
//     reader.read();
//
//     return 0;
// }
