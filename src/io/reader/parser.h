/**
 * @file parser.h
 * @brief Declares the common FEMaster input-deck parser lifecycle.
 *
 * `Parser` separates input processing according to semantic dependencies rather
 * than dense allocation requirements. Definition resources such as materials
 * and profiles are collected first so later part sections may resolve forward
 * references. Topology is then constructed in reusable parts and instances,
 * after which `Model::compile()` creates the dense solver representation.
 * Assembly-level sets/surfaces and enumeration-dependent fields are materialized
 * only after compilation; the final analysis pass executes loads, constraints
 * and load cases on compiled `ModelData`.
 *
 * The former identifier-counting stage is intentionally absent. Node, element
 * and surface identifiers remain sparse and part-local until compilation, so no
 * maximum-id capacities are required while parsing.
 *
 * @see ParserAbq
 * @see model::Model
 * @see model::Model::compile
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#pragma once

#include "../../core/types_num.h"
#include "../../loadcase/loadcase.h"
#include "../dsl/registry.h"
#include "../writer/writers.h"

#include <memory>
#include <string>

namespace fem {
namespace io { namespace dsl { class File; struct Line; } }
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
 * @brief Executes the common dependency-ordered input-deck workflow.
 *
 * The parser replays the complete deck in four dependency-ordered passes. Every
 * pass registers the same command grammar so nested scopes remain valid, but
 * `ActiveMode` selects which callbacks may mutate parser or model state:
 *
 * 1. collect global resources required by topology definitions,
 * 2. construct parts, instances and all part-local/default topology,
 * 3. materialize compiled assembly sets/surfaces plus dense fields,
 * 4. execute remaining loads, constraints and analysis commands.
 *
 * `Model::compile()` forms the explicit boundary between the topology and
 * assembly passes. Shell normals are completed at the end of the assembly pass
 * before any load case can access the compiled solver fields.
 *
 * `Parser` also owns the active load case while its consecutive commands are
 * interpreted, the result writers used during analysis and a separate registry
 * containing the complete grammar for documentation queries. Derived parsers
 * may specialize command activation for another deck dialect by overriding the
 * four pass-configuration functions.
 */
class Parser {
    // Persistent model and result-output state
    std::shared_ptr<model::Model> model_;
    io::writer::ResultWriters     writer_;

    // Complete command grammar used exclusively for documentation queries
    mutable io::dsl::Registry documentation_registry_;

    // Load case currently assembled by consecutive analysis-pass commands
    loadcase::LoadCase::Ptr active_loadcase_;
    int                     next_loadcase_id_ = 1;

public:
    // Construction
    Parser();
    virtual ~Parser();

    // Complete dependency-ordered deck evaluation and command documentation
    void run(const std::string&                   input_path,
             const std::string&                   output_path,
             const io::writer::WriterFileFormats& writer_formats = io::writer::WriterFileFormats());
    void document(const DocOptions& opts) const;

    // Current model and read-only documentation grammar
    const model::Model& model() const;
          model::Model& model();
    const io::dsl::Registry& registry() const;

    // Active load-case ownership. Activation assigns the next internal id and
    // common parser dependencies; completion executes and releases the case.
    void                begin_loadcase(loadcase::LoadCase::Ptr loadcase);
    void                end_loadcase();
    loadcase::LoadCase* active_loadcase();

protected:
    // Per-pass command activation. Derived deck readers may change which
    // callbacks execute while retaining the common four-pass lifecycle.
    virtual void configure_definition_pass(io::dsl::Registry& registry);
    virtual void configure_topology_pass  (io::dsl::Registry& registry);
    virtual void configure_assembly_pass  (io::dsl::Registry& registry);
    virtual void configure_analysis_pass  (io::dsl::Registry& registry);

private:
    // Individual complete-deck passes
    void run_definition_pass(const std::string& input_path);
    void run_topology_pass  (const std::string& input_path);
    void run_assembly_pass  (const std::string& input_path);
    void run_analysis_pass  (const std::string&                   input_path,
                             const std::string&                   output_path,
                             const io::writer::WriterFileFormats& writer_formats);

    // Complete FEMaster grammar and the persistent documentation view
    void register_commands(io::dsl::Registry& registry);
    void configure_documentation_registry();

};

} // namespace io::reader
} // namespace fem
