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

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        return 1;
    }

    namespace fs = std::filesystem;

    const int ncpus         = program.get<int>("--ncpus");
    std::string input_file  = program.get<std::string>("input_file");
    std::string output_file = program.get<std::string>("--output");
    const bool doc_mode     = program.get<bool>("--document");

    fem::global_config.max_threads = ncpus;

    // Single parser instance; ctor registers commands for doc mode immediately.
    fem::input_decks::Parser parser;

    if (doc_mode) {
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

        fem::input_decks::DocOptions opts;

        // Map format
        const auto fmt = program.get<std::string>("--doc-format");
        if      (fmt == "md")   opts.format = fem::input_decks::DocOptions::Format::Markdown;
        else if (fmt == "json") opts.format = fem::input_decks::DocOptions::Format::Json;
        else                    opts.format = fem::input_decks::DocOptions::Format::Text;

        // Map verbosity
        const auto verb = program.get<std::string>("--doc-verbosity");
        if      (verb == "index")   opts.verbosity = fem::input_decks::DocOptions::Verbosity::Index;
        else if (verb == "compact") opts.verbosity = fem::input_decks::DocOptions::Verbosity::Compact;
        else                        opts.verbosity = fem::input_decks::DocOptions::Verbosity::Full;

        // Wrap
        opts.wrap_width = program.get<int>("--doc-width");
        opts.no_wrap    = program.get<bool>("--doc-no-wrap");
        if (opts.no_wrap) opts.wrap_width = 0;

        // Regex (search)
        opts.regex = program.get<bool>("--doc-regex");

        // Action + payload
        if (program.get<bool>("--doc-list")) {
            opts.action = fem::input_decks::DocOptions::Action::List;
        } else if (program.is_used("--doc-show")) {
            opts.action = fem::input_decks::DocOptions::Action::Show;
            opts.cmd    = program.get<std::string>("--doc-show");
        } else if (program.is_used("--doc-tokens")) {
            opts.action = fem::input_decks::DocOptions::Action::Tokens;
            opts.cmd    = program.get<std::string>("--doc-tokens");
        } else if (program.is_used("--doc-variants")) {
            opts.action = fem::input_decks::DocOptions::Action::Variants;
            opts.cmd    = program.get<std::string>("--doc-variants");
        } else if (program.is_used("--doc-search")) {
            opts.action = fem::input_decks::DocOptions::Action::Search;
            opts.query  = program.get<std::string>("--doc-search");
        } else if (program.is_used("--doc-where-token")) {
            opts.action = fem::input_decks::DocOptions::Action::WhereToken;
            opts.query  = program.get<std::string>("--doc-where-token");
        } else if (program.get<bool>("--doc-all")) {
            opts.action = fem::input_decks::DocOptions::Action::All;
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

    // Regular solver mode
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

    fem::logging::info(true, "");
    fem::logging::info(true, "Input file : ", input_path.string());
    fem::logging::info(true, "Output file: ", output_file);
    fem::logging::info(true, "CPU(s)     : ", ncpus);
    fem::logging::info(true, "");

    try {
        parser.run(input_path.string(), output_file);
    } catch (const std::exception& e) {
        std::cerr << "Parsing failed: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
