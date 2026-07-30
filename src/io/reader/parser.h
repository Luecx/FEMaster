/**
 * @file parser.h
 * @brief Declares the staged FEMaster input-deck parser.
 *
 * `Parser` coordinates repeated passes over one input deck. The count pass
 * determines model capacities, the topology pass creates the discrete model,
 * the field pass creates fields that depend on element enumeration, and the
 * final data pass executes the remaining model and load-case commands.
 *
 * Command syntax and row parsing remain defined by the individual registration
 * helpers in `src/io/reader/commands`. Model-specific construction, including
 * shell-normal completion and equalisation, remains owned by `model::Model`.
 *
 * @see Parser
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
 * @struct DocOptions
 * @brief Configures generated documentation for registered DSL commands.
 *
 * The options select the documentation action, output representation and
 * verbosity. They do not affect normal model parsing or command execution.
 */
struct DocOptions {
    enum class Action { List, Show, Tokens, Variants, Search, WhereToken, All };
    enum class Format { Text, Markdown, Json };
    enum class Verbosity { Index, Compact, Full };

    // Requested documentation operation and output format
    Action    action    = Action::List;
    Format    format    = Format::Text;
    Verbosity verbosity = Verbosity::Full;

    // Optional command or search selection
    std::string cmd;
    std::string query;

    // Text rendering options
    int  wrap_width = 100;
    bool regex      = false;
    bool no_wrap    = false;
};

/**
 * @class Parser
 * @brief Executes the staged FEMaster input-deck workflow.
 *
 * The parser owns the active model, result writers, command registry and
 * load-case parsing state. Parsing is deliberately split into multiple passes
 * because field sizes depend on the element-node, integration-point and
 * material-point enumeration generated after topology construction.
 *
 * The field stage executes `*FIELD` and `*NORMAL` after element enumeration but
 * before `Model::build_shell_element_normals()`. Consequently, a selected
 * element-nodal field already exists when model-side shell normals are completed
 * and equalised. All other passes consume these commands without executing them
 * again.
 */
class Parser {
public:
    // Construction and destruction
    Parser();
    ~Parser();

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

private:
    /**
     * @struct CountData
     * @brief Stores the highest identifiers found during the allocation pass.
     */
    struct CountData {
        int highest_node_id    = -1;
        int highest_element_id = -1;
        int highest_surface_id = -1;
    };

    // Ordered parser stages
    CountData run_count_stage(const std::string& input_path);
    void run_topology_stage(const std::string& input_path);
    void run_field_stage(const std::string& input_path);
    void run_data_stage(const std::string&                   input_path,
                        const std::string&                   output_path,
                        const io::writer::WriterFileFormats& writer_formats);

    // Model allocation and command registration
    void allocate_model(const CountData& count);
    void register_count_commands(io::dsl::Registry& registry, CountData& count);
    void register_set_commands(io::dsl::Registry& registry);
    void register_topology_commands(io::dsl::Registry& registry);
    void register_analysis_commands(io::dsl::Registry& registry);
    void register_documentation_commands();

private:
    // Model and parser infrastructure
    std::shared_ptr<model::Model> m_model;
    io::writer::ResultWriters     m_writer;
    mutable io::dsl::Registry     m_registry;

    // Active load-case parsing state
    std::unique_ptr<loadcase::LoadCase> m_active_loadcase;
    std::string                         m_active_loadcase_type;
    int                                 m_next_loadcase_id = 1;
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
