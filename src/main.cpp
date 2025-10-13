#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <argparse/argparse.hpp>

#include "core/config.h"
#include "core/logging.h"
#include "input_decks/parser.h"

int main(int argc, char** argv) {
    argparse::ArgumentParser program("FEM Solver");

    // Optional input file (only required if not using --document)
    program.add_argument("input_file")
        .default_value(std::string{})   // allow running without it
        .help("Path to the input file (.inp default extension). Optional when using --document.");

    program.add_argument("--ncpus")
        .default_value(1)
        .scan<'i', int>()
        .help("Number of CPUs to use (default: 1)");

    program.add_argument("--output")
        .default_value(std::string{})
        .help("Override output filename (default: input with .res extension)");

    program.add_argument("--document")
        .default_value(std::string{})
        .help("Print DSL documentation for KEYWORD (use ALL for full list)");

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        return 1;
    }

    namespace fs = std::filesystem;

    int ncpus = program.get<int>("ncpus");
    std::string input_file  = program.get<std::string>("input_file");
    std::string output_file = program.get<std::string>("output");
    std::string doc_keyword = program.get<std::string>("document");

    // Require input_file if not using --document
    if (input_file.empty() && doc_keyword.empty()) {
        std::cerr << "Error: Missing input file. "
                     "You must provide one unless --document is specified."
                  << std::endl;
        return 1;
    }

    // Handle documentation-only mode
    if (!doc_keyword.empty()) {
        fem::logging::info(true, "");
        fem::logging::info(true, "FEMaster Documentation Mode");
        fem::logging::info(true, "Keyword: ", doc_keyword);
        fem::logging::info(true, "");

        try {
            fem::input_decks::Parser parser("", "");  // no real files
            parser.prepare_for_documentation();
            parser.register_all_commands();

            if (doc_keyword == "ALL")
                parser.registry().print_help();
            else
                parser.registry().print_help(doc_keyword);

        } catch (const std::exception& e) {
            std::cerr << "Documentation generation failed: " << e.what() << std::endl;
            return 1;
        }
        return 0;
    }

    // Regular solver mode
    if (input_file.find(".inp") == std::string::npos)
        input_file += ".inp";

    fs::path input_path = fs::path(input_file);
    if (!fs::exists(input_path)) {
        std::cerr << "Error: Input file '" << input_path.string() << "' does not exist." << std::endl;
        return 1;
    }

    if (output_file.empty()) {
        fs::path out = input_path;
        out.replace_extension(".res");
        output_file = out.string();
    }

    fem::logging::info(true, "");
    fem::logging::info(true, "Input file : ", input_path.string());
    fem::logging::info(true, "Output file: ", output_file);
    fem::logging::info(true, "CPU(s)     : ", ncpus);
    fem::logging::info(true, "");

    fem::global_config.max_threads = ncpus;

    try {
        fem::input_decks::Parser parser(input_path.string(), output_file);
        parser.preprocess();
        parser.register_all_commands();
        parser.parse();
    } catch (const std::exception& e) {
        std::cerr << "Parsing failed: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
