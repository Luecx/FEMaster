/**
 * @file parser.h
 * @brief Declares the common staged FEMaster input-deck parser.
 *
 * `Parser` coordinates repeated passes over one input deck. The count pass
 * determines model capacities, the topology pass creates the discrete model,
 * the field pass creates fields that depend on element enumeration, and the
 * final data pass executes the remaining model and load-case commands.
 *
 * The native FEMaster reader provides the default command registration for all
 * four stages. Syntax specializations may override only the stage-specific
 * command registration and activation while reusing model allocation, section
 * assignment, element-local enumeration, shell-normal construction and result
 * writer setup.
 *
 * Command syntax and row parsing remain defined by the individual registration
 * helpers in `src/io/reader/commands` and format-specific command directories.
 * Model construction and finite-element semantics remain owned by `model::Model`
 * and the corresponding element, section, material and load-case abstractions.
 *
 * @see ParserAbq
 * @see model::Model
 *
 * @author Finn Eggers
 * @date 30.07.2026
 */

#pragma once

#include "../../core/types_num.h"
#include "../dsl/registry.h"
#include "../writer/writers.h"

#include <memory>
#include <string>

namespace fem {
namespace loadcase { struct LoadCase; }
namespace io { namespace dsl { class File; struct Line; } }
namespace model { struct Model; }

namespace io::reader {

/**
 * @brief Configures generated documentation for registered DSL commands.
 *
 * The options select one documentation operation together with its optional
 * command/search payload, text representation and wrapping behavior. They are
 * used only by documentation mode and do not modify normal model parsing or
 * solver execution.
 */
struct DocOptions {
    // Documentation modes
    enum class Action { List, Show, Tokens, Variants, Search, WhereToken, All };
    enum class Format { Text, Markdown, Json };
    enum class Verbosity { Index, Compact, Full };

    // Requested documentation operation
    Action action = Action::List;

    // Optional command or search selection
    std::string cmd;
    std::string query;

    // Output representation and text rendering
    Format    format     = Format::Text;
    Verbosity verbosity  = Verbosity::Full;
    int       wrap_width = 100;
    bool      regex      = false;
    bool      no_wrap    = false;
};

/**
 * @brief Executes the common dependency-ordered input-deck workflow.
 *
 * The parser owns the active model, result writers, documentation registry and
 * load-case parsing state. Parsing is split into four repeated passes because
 * model capacities must be known before allocation and generic field dimensions
 * depend on element-local enumeration generated after topology construction.
 *
 * The fixed execution order is:
 *
 * 1. count identifiers required for model allocation,
 * 2. construct topology and assign/enumerate element-local data,
 * 3. create enumeration-dependent fields and select shell-normal data,
 * 4. execute all remaining model and analysis commands.
 *
 * Derived readers preserve this sequence and override only the registry
 * configuration of each stage. This keeps format-specific syntax separate from
 * FEMaster model construction and solver behavior.
 */
class Parser {
protected:
    /**
     * @brief Allocation information collected during the count pass.
     *
     * The parser allocates dense FEMaster repositories from the highest external
     * identifiers rather than the number of encountered records. Unused entity
     * classes remain at `-1`, producing zero-capacity domains after the common
     * `+1` allocation conversion.
     *
     * Derived readers update the same counters from their format-specific count
     * registrations so all syntax variants share identical model allocation.
     */
    struct CountData {
        int highest_node_id    = -1;
        int highest_element_id = -1;
        int highest_surface_id = -1;
    };

private:
    // Model and parser infrastructure shared by every syntax specialization
    std::shared_ptr<model::Model> m_model;
    io::writer::ResultWriters     m_writer;
    mutable io::dsl::Registry     m_registry;

    // Active load-case parsing state used by native command callbacks
    std::unique_ptr<loadcase::LoadCase> m_active_loadcase;
    std::string                         m_active_loadcase_type;
    int                                 m_next_loadcase_id = 1;

public:
    // Construction and destruction
    Parser();
    virtual ~Parser();

    // High-level parsing and documentation operations
    void run(const std::string&                   input_path,
             const std::string&                   output_path,
             const io::writer::WriterFileFormats& writer_formats = io::writer::WriterFileFormats());
    void document(const DocOptions& opts) const;

    // Parsed model, writer and documentation registry access
    model::Model& model();
    const model::Model& model() const;
    io::writer::ResultWriters& writer();
    const io::writer::ResultWriters& writer() const;
    io::dsl::Registry& registry();
    const io::dsl::Registry& registry() const;

    // Load-case state used by command registration callbacks
    int next_loadcase_id();
    void set_active_loadcase(std::unique_ptr<loadcase::LoadCase> lc, std::string type);
    loadcase::LoadCase* active_loadcase();
    const loadcase::LoadCase* active_loadcase() const;
    template<class T> T* active_loadcase_as();
    template<class T> const T* active_loadcase_as() const;
    void clear_active_loadcase();
    const std::string& active_loadcase_type() const;

protected:
    // Format-specific command registration and activation for the four common
    // passes. Derived readers may change syntax handling but not pass ordering.
    virtual void configure_count_stage   (io::dsl::Registry& registry, CountData& count);
    virtual void configure_topology_stage(io::dsl::Registry& registry);
    virtual void configure_field_stage   (io::dsl::Registry& registry);
    virtual void configure_data_stage    (io::dsl::Registry& registry);

private:
    // Ordered parser stages and model allocation
    CountData run_count_stage(const std::string& input_path);
    void run_topology_stage(const std::string& input_path);
    void run_field_stage(const std::string& input_path);
    void run_data_stage(const std::string&                   input_path,
                        const std::string&                   output_path,
                        const io::writer::WriterFileFormats& writer_formats);
    void allocate_model(const CountData& count);

    // Native FEMaster command registration used by the default stage
    // configuration and documentation mode.
    void register_count_commands(io::dsl::Registry& registry, CountData& count);
    void register_set_commands(io::dsl::Registry& registry);
    void register_topology_commands(io::dsl::Registry& registry);
    void register_analysis_commands(io::dsl::Registry& registry);
    void register_documentation_commands();
};

// Typed access to the active load case

template<class T>
inline T* Parser::active_loadcase_as() {
    return dynamic_cast<T*>(active_loadcase());
}

template<class T>
inline const T* Parser::active_loadcase_as() const {
    return dynamic_cast<const T*>(active_loadcase());
}

} // namespace io::reader
} // namespace fem
