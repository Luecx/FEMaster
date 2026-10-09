/**
 * @file main.cpp
 * @brief Defines the FEMaster command-line entry point and unified input parsing.
 *
 * The executable parses command-line options, configures global runtime settings
 * and passes both native FEMaster and supported Abaqus syntax through one parser.
 * The --format switch remains available as a backwards-compatible CLI alias.
 *
 * Documentation mode uses the unified FEMaster/Abaqus command registry and
 * is handled before regular input-file validation and parser dispatch.
 *
 * @see fem::io::reader::Parser
 *
 * @author Finn Eggers
 * @date 17.08.2026
 */

#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <argparse/argparse.hpp>

#include "core/config.h"
#include "core/timer.h"
#include "core/logging.h"
#include "core/version.h"
#include "io/reader/parser.h"

/**
 * Runs FEMaster in documentation or solver mode.
 *
 * Command-line parsing is completed before any model reader is constructed.
 * Documentation requests use the unified FEMaster parser registry, while regular
 * solver runs validate the input/output paths and invoke the unified parser.
 * The --format alias, writer settings and thread limits are independent.
 *
 * @param argc Number of command-line arguments.
 * @param argv Command-line argument values.
 * @return Zero on successful completion and one after argument, documentation,
 *         input-file or parser errors.
 */
int main(int argc, char** argv) {
    auto timer = fem::Timer();
    timer.start();

    // Use the same solver version for --version and the startup banner.
    const std::string version = std::to_string(VERSION_MAJOR) + "." +
                                std::to_string(VERSION_MINOR) + "." +
                                std::to_string(VERSION_PATCH);
    argparse::ArgumentParser program("FEM Solver", version);

    // Optional input file (only required if not in doc mode)
    program.add_argument("input_file")
        .default_value(std::string{})
        .help("Path to the input file (.inp default extension). Optional in documentation mode.");

    program.add_argument("--ncpus")
        .default_value(1)
        .scan<'i', int>()
        .help("Number of CPUs to use (default: 1)");

    program.add_argument("--output")
        .default_value(std::string{})
        .help("Override output filename (default: input with .res extension)");

    program.add_argument("--format")
        .default_value(std::string{"femaster"})
        .choices("femaster", "abaqus")
        .help("Input format alias (femaster or abaqus); both use the unified parser.");

    program.add_argument("--output-format")
        .nargs(argparse::nargs_pattern::at_least_one)
        .choices("res", "frd", "femr")
        .help("Result formats to write (default: res frd).");

    // ---- Documentation mode flags (flat, no nested parser) ----
    program.add_argument("--document")
        .flag()
        .help("Enter documentation mode (ignores input/output run). Choose one --doc-* action.");

    // Actions (exactly one when --document is set)
    program.add_argument("--doc-list").flag().help("List all commands (index).");
    program.add_argument("--doc-show").help("Show full docs for one command.");
    program.add_argument("--doc-tokens").help("List tokens for one command.");
    program.add_argument("--doc-variants").help("List variants for one command.");
    program.add_argument("--doc-search").help("Search across command names and descriptions.");
    program.add_argument("--doc-where-token").help("Find commands containing a token (by name/desc).");
    program.add_argument("--doc-all").flag().help("Print full documentation for all commands.");

    // Formatting knobs
    program.add_argument("--doc-format")
        .default_value(std::string{"text"})
        .choices("text","md","json")
        .help("Output format for documentation (text supported now).");

    program.add_argument("--doc-verbosity")
        .default_value(std::string{"full"})
        .choices("index","compact","full")
        .help("Verbosity for --doc-show.");

    program.add_argument("--doc-width")
        .scan<'i', int>()
        .default_value(100)
        .help("Wrap width for docs (0 = no wrap).");

    program.add_argument("--doc-no-wrap")
        .flag()
        .help("Disable wrapping (same as --doc-width 0).");

    program.add_argument("--doc-regex")
        .flag()
        .help("Treat --doc-search as regex.");

    // Parse and validate the command-line syntax before constructing a reader
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        return 1;
    }

    namespace fs = std::filesystem;

    const int         ncpus    = program.get<int>("--ncpus");
    const std::string format   = program.get<std::string>("--format");
    const bool        doc_mode = program.get<bool>("--document");
    std::string       input_file  = program.get<std::string>("input_file");
    std::string       output_file = program.get<std::string>("--output");

    fem::io::writer::WriterFileFormats writer_formats;
    if (program.is_used("--output-format")) {
        writer_formats = {};
        writer_formats.res = false;
        writer_formats.frd = false;
        for (const auto& output_format : program.get<std::vector<std::string>>("--output-format")) {
            writer_formats.res  = writer_formats.res  || output_format == "res";
            writer_formats.frd  = writer_formats.frd  || output_format == "frd";
            writer_formats.femr = writer_formats.femr || output_format == "femr";
        }
    }

    fem::global_config.max_threads = ncpus;

    // Handle native command documentation independently of solver input syntax
    if (doc_mode) {
        fem::io::reader::Parser parser;

        // Enforce exactly one action
        int actions = 0;
        actions += program.get<bool>("--doc-list") ? 1 : 0;
        actions += program.is_used("--doc-show") ? 1 : 0;
        actions += program.is_used("--doc-tokens") ? 1 : 0;
        actions += program.is_used("--doc-variants") ? 1 : 0;
        actions += program.is_used("--doc-search") ? 1 : 0;
        actions += program.is_used("--doc-where-token") ? 1 : 0;
        actions += program.get<bool>("--doc-all") ? 1 : 0;

        if (actions != 1) {
            std::cerr << "In doc mode, choose exactly ONE action among: "
                         "--doc-list | --doc-show CMD | --doc-tokens CMD | --doc-variants CMD | "
                         "--doc-search TEXT | --doc-where-token TOKEN | --doc-all\n";
            return 1;
        }

        fem::io::reader::DocOptions opts;

        // Map format
        const auto fmt = program.get<std::string>("--doc-format");
        if      (fmt == "md")   opts.format = fem::io::reader::DocOptions::Format::Markdown;
        else if (fmt == "json") opts.format = fem::io::reader::DocOptions::Format::Json;
        else                    opts.format = fem::io::reader::DocOptions::Format::Text;

        // Map verbosity
        const auto verb = program.get<std::string>("--doc-verbosity");
        if      (verb == "index")   opts.verbosity = fem::io::reader::DocOptions::Verbosity::Index;
        else if (verb == "compact") opts.verbosity = fem::io::reader::DocOptions::Verbosity::Compact;
        else                        opts.verbosity = fem::io::reader::DocOptions::Verbosity::Full;

        // Wrap
        opts.wrap_width = program.get<int>("--doc-width");
        opts.no_wrap    = program.get<bool>("--doc-no-wrap");
        if (opts.no_wrap) opts.wrap_width = 0;

        // Regex (search)
        opts.regex = program.get<bool>("--doc-regex");

        // Action + payload
        if (program.get<bool>("--doc-list")) {
            opts.action = fem::io::reader::DocOptions::Action::List;
        } else if (program.is_used("--doc-show")) {
            opts.action = fem::io::reader::DocOptions::Action::Show;
            opts.cmd    = program.get<std::string>("--doc-show");
        } else if (program.is_used("--doc-tokens")) {
            opts.action = fem::io::reader::DocOptions::Action::Tokens;
            opts.cmd    = program.get<std::string>("--doc-tokens");
        } else if (program.is_used("--doc-variants")) {
            opts.action = fem::io::reader::DocOptions::Action::Variants;
            opts.cmd    = program.get<std::string>("--doc-variants");
        } else if (program.is_used("--doc-search")) {
            opts.action = fem::io::reader::DocOptions::Action::Search;
            opts.query  = program.get<std::string>("--doc-search");
        } else if (program.is_used("--doc-where-token")) {
            opts.action = fem::io::reader::DocOptions::Action::WhereToken;
            opts.query  = program.get<std::string>("--doc-where-token");
        } else if (program.get<bool>("--doc-all")) {
            opts.action = fem::io::reader::DocOptions::Action::All;
        }

        try {
            fem::logging::info(true, "");
            fem::logging::info(true, "FEMaster Documentation Mode");
            fem::logging::info(true, "");
            parser.document(opts);
        } catch (const std::exception& e) {
            std::cerr << "Documentation generation failed: " << e.what() << std::endl;
            return 1;
        }
        return 0;
    }

    // Validate and normalize regular solver input/output paths
    if (input_file.empty()) {
        std::cerr << "Error: Missing input file. You must provide one unless --document is specified.\n";
        return 1;
    }
    if (input_file.find(".inp") == std::string::npos)
        input_file += ".inp";

    fs::path input_path = fs::path(input_file);
    if (!fs::exists(input_path)) {
        std::cerr << "Error: Input file '" << input_path.string() << "' does not exist.\n";
        return 1;
    }

    if (output_file.empty()) {
        fs::path out = input_path;
        out.replace_extension(".res");
        output_file = out.string();
    }

    // Report the resolved run configuration before parsing the model
    fem::logging::info(true, "");
    fem::logging::info(true, "Input file : ", input_path.string());
    fem::logging::info(true, "Output file: ", output_file);
    fem::logging::info(true, "Format     : ", format);
    fem::logging::info(true, "CPU(s)     : ", ncpus);
    fem::logging::info(true, "Write .res : ", writer_formats.res ? "yes" : "no");
    fem::logging::info(true, "Write .frd : ", writer_formats.frd ? "yes" : "no");
    fem::logging::info(true, "Write .femr: ", writer_formats.femr ? "yes" : "no");
    fem::logging::info(true, "");

    // Both format aliases use the same native/Abaqus-compatible grammar.
    try {
        fem::io::reader::Parser parser;
        parser.run(input_path.string(), output_file, writer_formats);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    // Report elapsed wall time in h, min, sec and ms
    timer.stop();
    auto t  = timer.elapsed();
    auto ms = t % 1000;
    auto s  = t / 1000 % 60;
    auto m  = t / 60000 % 60;
    auto h  = t / 3600000;

    fem::logging::info(
        true,
        "Process finished in ",
        h ? std::to_string(h) + "h " : "",
        m ? std::to_string(m) + "min " : "",
        s ? std::to_string(s) + "s " : "",
        ms, "ms. (", t, " ms total)"
    );

    return 0;
}
