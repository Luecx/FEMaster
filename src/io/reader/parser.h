/**
 * @file parser.h
 * @brief Declares the common FEMaster input-deck parser lifecycle.
 *
 * The reader now separates syntax from semantic model construction. A complete
 * command grammar is registered once, the source deck is parsed once into a
 * reusable `dsl::Deck`, and semantic commands are then executed explicitly in
 * dependency order. `Model::compile()` remains the visible one-way boundary
 * between Part/Instance construction and assembly-level materialization.
 *
 * This removes parser stages and command activation modes entirely from the
 * reader lifecycle. Scope-dependent commands such as `NSET`, `ELSET`, `SURFACE`
 * and point-element properties are parsed once and later selected by their
 * stored parent occurrence when the required model state exists.
 *
 * @see ParserAbq
 * @see io::dsl::Deck
 * @see io::dsl::DeckParser
 * @see model::Model::compile
 *
 * @author Finn Eggers
 * @date 26.08.2026
 */

#pragma once

#include "../../core/types_num.h"
#include "../../loadcase/loadcase.h"
#include "../dsl/deck.h"
#include "../dsl/registry.h"
#include "../writer/writers.h"

#include <memory>
#include <string>

namespace fem {
namespace model { struct Model; }

namespace io::reader {

/**
 * @brief Selects and formats command-language documentation output.
 *
 * The action determines which registry view is printed and which command or
 * search query is used. Text is currently the implemented output format;
 * Markdown, JSON and wrapping controls are accepted by the option contract but
 * are not yet applied by the current text renderer.
 */
struct DocOptions {
    // Documentation operation and output representation
    enum class Action { List, Show, Tokens, Variants, Search, WhereToken, All };
    enum class Format { Text, Markdown, Json };
    enum class Verbosity { Index, Compact, Full };

    // Selected operation and optional command/search arguments
    Action action = Action::List;

    std::string cmd;
    std::string query;

    // Output formatting controls
    Format    format     = Format::Text;
    Verbosity verbosity  = Verbosity::Full;
    int       wrap_width = 100;
    bool      regex      = false;
    bool      no_wrap    = false;
};

/**
 * @brief Parses one deck once and executes its semantic commands explicitly.
 *
 * `Parser::run()` performs only the common lifecycle: reset state, register the
 * complete grammar, parse the file into a reusable tree and delegate semantic
 * construction to `process_deck()`. The base FEMaster reader and specialized
 * Abaqus reader provide their own explicit command order while sharing the same
 * parsed-deck representation and load-case ownership.
 *
 * `Model::compile()` is called from semantic processing exactly where the model
 * crosses from sparse Part/Instance topology to dense assembly storage. No
 * command activation state is stored in the DSL registry.
 */
class Parser {
    // Persistent model and result-output state
    std::shared_ptr<model::Model> model_;
    io::writer::ResultWriters     writer_;

    // Complete command grammar used exclusively for documentation queries
    mutable io::dsl::Registry documentation_registry_;

    // Load case currently assembled by consecutive analysis commands
    loadcase::LoadCase::Ptr active_loadcase_;
    int                     next_loadcase_id_ = 1;

public:
    // Construction
    Parser();
    virtual ~Parser();

    // Parse once, process explicitly and expose command documentation
    void run(const std::string&                   input_path,
             const std::string&                   output_path,
             const io::writer::WriterFileFormats& writer_formats = io::writer::WriterFileFormats());
    void document(const DocOptions& opts) const;

    // Current model and read-only documentation grammar
    const model::Model& model() const;
          model::Model& model();
    const io::dsl::Registry& registry() const;

    // Active load-case ownership used by LOADCASE/STEP command callbacks
    void                begin_loadcase(loadcase::LoadCase::Ptr loadcase);
    void                end_loadcase();
    loadcase::LoadCase* active_loadcase();

protected:
    // Dialect-specific grammar and explicit semantic processing order
    virtual void register_commands(io::dsl::Registry& registry);
    virtual void process_deck(const io::dsl::Deck&                  deck,
                              const std::string&                    input_path,
                              const std::string&                    output_path,
                              const io::writer::WriterFileFormats& writer_formats);

    // Result writers are shared by the native and Abaqus processing paths
    void initialize_writers(const std::string&                    input_path,
                            const std::string&                    output_path,
                            const io::writer::WriterFileFormats& writer_formats);
    void close_writers();

    // Rebuild documentation after the complete dialect grammar is known
    void configure_documentation_registry();
};

} // namespace io::reader
} // namespace fem
