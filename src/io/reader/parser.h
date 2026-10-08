/**
 * @file parser.h
 * @brief Declares the unified FEMaster and Abaqus input-deck parser lifecycle.
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
 * @see io::dsl::Deck
 * @see io::dsl::DeckParser
 * @see model::Model::compile
 *
 * @author Finn Eggers
 * @date 26.08.2026
 */

#pragma once

#include "../../bc/condition_manager.h"
#include "../../core/types_num.h"
#include "../../loadcase/loadcase.h"
#include "../dsl/deck.h"
#include "../dsl/registry.h"
#include "../writer/writers.h"

#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <utility>

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
 * complete grammar, parse the file into a reusable tree and execute semantic
 * construction in `process_deck()`. Both supported input dialects share this
 * grammar, command order and load-case ownership.
 *
 * `Model::compile()` is called from semantic processing exactly where the model
 * crosses from sparse Part/Instance topology to dense assembly storage. No
 * command activation state is stored in the DSL registry.
 *
 * Input-condition history is resolved here, not in ConditionManager. Each
 * history family has two identifier maps tracking inherited and current names
 * to one or more physical Condition pointers. Input identifiers are never
 * stored on Condition objects. One CLOAD may create multiple physical pointers
 * when nodes use different TRANSFORM coordinate systems.
 * The parser resolves MOD using the original source and DOF. Inherited and
 * current-step identifiers are kept separately; the manager owns only active
 * physical Condition pointers.
 *
 * Named collector definitions are independent of direct input history.
 * Selecting a named collector records the shared physical definitions for
 * the current analysis. When that analysis executes, only pointers that were
 * not already active are inserted temporarily and removed afterward.
 * This preserves direct conditions that happen to share a collector pointer.
 */
/**
 * Step syntax is a scope; its procedure constructs the actual load case.
 * These values are step-level controls, never separate condition storage.
 */
struct StepState {
    bool        step_active    = false;
    int         max_increments = 100;
    bool        nlgeom         = false;
    bool        perturbation   = false;
    Precision   step_period    = Precision(1);
};

class Parser {
    // Persistent model and result-output state
    std::shared_ptr<model::Model> model_;
    io::writer::ResultWriters     writer_;

    // Complete command grammar used exclusively for documentation queries
    mutable io::dsl::Registry documentation_registry_;

    // Load case currently assembled by consecutive analysis commands
    loadcase::LoadCase::Ptr active_loadcase_;
    int                     next_loadcase_id_ = 1;
    StepState               step_state_;

    // Reader-owned input identity, separate from active physical conditions.
    // One identifier may resolve to multiple TRANSFORM fragments.
    using IdentifierMap = std::unordered_map<std::string, std::vector<bc::Condition::Ptr>>;
    std::array<IdentifierMap, bc::N_CONDITION_FAMILIES> inherited_conditions_;
    std::array<IdentifierMap, bc::N_CONDITION_FAMILIES> current_conditions_;

    // Shared collector definitions requested by the current analysis. These
    // are inserted into the manager only while that analysis is running.
    std::vector<std::pair<bc::ConditionFamily, bc::Condition::Ptr>> selected_conditions_;

public:
    // we allow *TRANSFORM to be parsed but we don't really transform the DOFs. Instead,
    // we apply the transformation when creating supports and loads. That way, we can handle
    // two loads at the same node in different coordinate systems. We need to track all nodes
    // that have been transformed so that, once we encounter a node, we can create a separate
    // load with the correct transformation applied. This maps node ID -> coordinate system name.
    std::unordered_map<ID, std::string> node_transforms;

public:
    // Construction
    Parser();
    ~Parser();

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
    StepState& step_state();
    const StepState& step_state() const;

    // Physical conditions remain model-owned; input identities remain here.
    bc::ConditionManager&       conditions();
    const bc::ConditionManager& conditions() const;

    // MOD replaces one complete logical source/DOF group. Compatible physical
    // fragments retain their start values; changed bases keep both histories.
    void add_condition    (bc::ConditionFamily family, bc::Condition::Ptr condition,
                            const std::string& identifier = {});
    void modify_conditions(bc::ConditionFamily family, const std::string& identifier,
                            std::vector<bc::Condition::Ptr> replacements);
    void clear_conditions (bc::ConditionFamily family);

    // Select a reusable collector definition for the current analysis only.
    // Repeated selection of the same pointer and family is idempotent.
    void select_collector_condition(bc::ConditionFamily family, bc::Condition::Ptr condition);

    // Activate an explicitly named definition target; never create implicit collectors
    std::string activate_load_collector(const std::string& name);

protected:
    // Unified grammar and explicit semantic processing order
    void register_commands(io::dsl::Registry& registry);
    void process_deck(const io::dsl::Deck&                  deck,
                              const std::string&                    input_path,
                              const std::string&                    output_path,
                              const io::writer::WriterFileFormats& writer_formats);

    // Result writers are shared by all supported input syntax
    void initialize_writers(const std::string&                    input_path,
                            const std::string&                    output_path,
                            const io::writer::WriterFileFormats& writer_formats);
    void close_writers();

    // Rebuild documentation after the complete dialect grammar is known
    void configure_documentation_registry();
};

} // namespace io::reader
} // namespace fem
