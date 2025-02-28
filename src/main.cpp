
#include "reader/reader.h"
#include "core/logging.h"

#include <argparse/argparse.hpp>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <filesystem>


#include "math/interpolate.h"

int main(int argc, char* argv[]) {
    // Set up argparse
    argparse::ArgumentParser program("FEM Solver");

    program.add_argument("input_file")
        .help("Path to the input file (.inp is optional).");

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

    // Ensure the input file has ".inp" extension
    if (input_file.find(".inp") == std::string::npos) {
        input_file += ".inp";
    }

    // Check if input file exists
    if (!std::filesystem::exists(input_file)) {
        std::cerr << "Error: Input file '" << input_file << "' does not exist." << std::endl;
        return 1;
    }

    // Create the output file by replacing ".inp" with ".res"
    std::string output_file = input_file;
    size_t pos = output_file.find(".inp");
    if (pos != std::string::npos) {
        output_file.replace(pos, 4, ".res");
    }

    // Logging input and output file information
    fem::logging::info(true, "");
    fem::logging::info(true, "Input file : ", input_file);
    fem::logging::info(true, "Output file: ", output_file);
    fem::logging::info(true, "CPU(s)     : ", ncpus);
    fem::logging::info(true, "");

    // Store number of CPUs in config
    fem::global_config.max_threads = ncpus;

    // Read the input file using the reader
    fem::reader::Reader reader{input_file, output_file};
    reader.read();

    return 0;
}

// #include "material/elasticity.h"
// #include "material/orthotropic_elasticity.h"
// #include "material/material.h"
//
// #include "cos/rectangular_system.h"
// #include "cos/coordinate_system.h"
//
// #include <iostream>
// #include <Eigen/Dense>
// #include <cmath>
//
// using namespace fem;
//
// Precision epsilon = 1e-5;  // Perturbation for numerical differentiation
//
// // Function to compute the transformed material matrix numerically
// StaticMatrix<6,6> numerical_transformed_derivative(fem::material::Elasticity& elasticity, Precision angle_1, Precision angle_2, Precision angle_3, int angle_index) {
//     StaticMatrix<3,3> R_plus, R_minus;
//
//     if (angle_index == 0) {
//         R_plus = cos::RectangularSystem::euler(angle_1 + epsilon, angle_2, angle_3).get_axes(Vec3(0,0,0));
//         R_minus = cos::RectangularSystem::euler(angle_1 - epsilon, angle_2, angle_3).get_axes(Vec3(0,0,0));
//     } else if (angle_index == 1) {
//         R_plus = cos::RectangularSystem::euler(angle_1, angle_2 + epsilon, angle_3).get_axes(Vec3(0,0,0));
//         R_minus = cos::RectangularSystem::euler(angle_1, angle_2 - epsilon, angle_3).get_axes(Vec3(0,0,0));
//     } else if (angle_index == 2) {
//         R_plus = cos::RectangularSystem::euler(angle_1, angle_2, angle_3 + epsilon).get_axes(Vec3(0,0,0));
//         R_minus = cos::RectangularSystem::euler(angle_1, angle_2, angle_3 - epsilon).get_axes(Vec3(0,0,0));
//     }
//
//     StaticMatrix<6,6> transformed_plus = elasticity.get_transformed<3>(R_plus);
//     StaticMatrix<6,6> transformed_minus = elasticity.get_transformed<3>(R_minus);
//
//     return (transformed_plus - transformed_minus) / (2 * epsilon);
// }
//
// int main() {
//     Precision E = 210;
//     Precision v = 0.3;
//     Precision G = E / (2 * (1 + v));
//
//     fem::material::OrthotropicElasticity ortho_elasticity{E, E, E, G, G, G, v, v, v};
//
//     Precision angle_1 = 0.1;
//     Precision angle_2 = 0.2;
//     Precision angle_3 = 0.3;
//
//     // // Rotation matrix and its derivatives
//     // auto transform = cos::RectangularSystem::euler(angle_1, angle_2, angle_3).get_axes(Vec3(0,0,0));
//     // auto transform_p1 = cos::RectangularSystem::euler(angle_1 + epsilon, angle_2, angle_3).get_axes(Vec3(0,0,0));
//     // auto transform_m1 = cos::RectangularSystem::euler(angle_1 - epsilon, angle_2, angle_3).get_axes(Vec3(0,0,0));
//     //
//     // auto derivative_1 = cos::RectangularSystem::derivative_rot_x(angle_1, angle_2, angle_3);
//     //
//     // std::cout << transform << std::endl;
//     // std::cout << std::endl;
//     // std::cout << derivative_1 << std::endl;
//     // std::cout << std::endl;
//     // std::cout << (transform_p1 - transform_m1) / (2 * epsilon) << std::endl;
//     // std::cout << std::endl;
//     //
//     // auto mat_transform_der_1 = ortho_elasticity.transformation_der(transform, derivative_1);
//     // auto mat_transform_der_1_num =  (ortho_elasticity.transformation(transform_p1) - ortho_elasticity.transformation(transform_m1)) / (2 * epsilon);
//     //
//     // std::cout << "Analytical Derivative w.r.t. angle_1:\n" << mat_transform_der_1 << std::endl;
//     // std::cout << "Numerical Derivative w.r.t. angle_1:\n" << mat_transform_der_1_num << std::endl;
//
//     // auto derivative_2 = cos::RectangularSystem::derivative_rot_y(angle_1, angle_2, angle_3);
//     // auto derivative_3 = cos::RectangularSystem::derivative_rot_z(angle_1, angle_2, angle_3);
//
//     // // Compute analytical derivatives
//     // StaticMatrix<6,6> analytical_derivative_1 = ortho_elasticity.get_transformed_derivative<3>(transform, derivative_1);
//     // // StaticMatrix<6,6> analytical_derivative_2 = ortho_elasticity.get_transformed_derivative<3>(transform, derivative_2);
//     // // StaticMatrix<6,6> analytical_derivative_3 = ortho_elasticity.get_transformed_derivative<3>(transform, derivative_3);
//     //
//     // // Compute numerical derivatives
//     // StaticMatrix<6,6> numerical_derivative_1 = numerical_transformed_derivative(ortho_elasticity, angle_1, angle_2, angle_3, 0);
//     // // StaticMatrix<6,6> numerical_derivative_2 = numerical_transformed_derivative(ortho_elasticity, angle_1, angle_2, angle_3, 1);
//     // // StaticMatrix<6,6> numerical_derivative_3 = numerical_transformed_derivative(ortho_elasticity, angle_1, angle_2, angle_3, 2);
//     //
//     // Output results
//     // std::cout << "Analytical Derivative w.r.t. angle_1:\n" << analytical_derivative_1 << std::endl;
//     // std::cout << "Numerical Derivative w.r.t. angle_1:\n" << numerical_derivative_1 << std::endl;
//     // std::cout << "Analytical Derivative w.r.t. angle_2:\n" << analytical_derivative_2 << std::endl;
//     // std::cout << "Numerical Derivative w.r.t. angle_2:\n" << numerical_derivative_2 << std::endl;
//     // std::cout << "Analytical Derivative w.r.t. angle_3:\n" << analytical_derivative_3 << std::endl;
//     // std::cout << "Numerical Derivative w.r.t. angle_3:\n" << numerical_derivative_3 << std::endl;
//
//     return 0;
// }
