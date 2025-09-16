
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
