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
 * only after compilation; the final data pass executes loads, constraints and
 * analyses on compiled `ModelData`.
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
#include "../dsl/registry.h"
#include "../writer/writers.h"

#include <memory>
#include <string>

namespace fem {
namespace loadcase { struct LoadCase; }
namespace io { namespace dsl { class File; struct Line; } }
namespace model { struct Model; }

namespace io::reader {

struct DocOptions {
    enum class Action { List, Show, Tokens, Variants, Search, WhereToken, All };
    enum class Format { Text, Markdown, Json };
    enum class Verbosity { Index, Compact, Full };

    Action action = Action::List;

    std::string cmd;
    std::string query;

    Format    format     = Format::Text;
    Verbosity verbosity  = Verbosity::Full;
    int       wrap_width = 100;
    bool      regex      = false;
    bool      no_wrap    = false;
};

/**
 * @brief Executes the common dependency-ordered input-deck workflow.
 *
 * Model construction is split into four semantic passes:
 *
 * 1. collect global resources required by topology definitions,
 * 2. construct parts, instances and all part-local/default topology,
 * 3. compile and materialize assembly sets/surfaces plus dense fields,
 * 4. execute remaining loads, constraints and analysis commands.
 *
 * The historical `field` method names are retained internally for the third
 * pass, but that pass is now the complete post-compile materialization stage.
 */
class Parser {
private:
    std::shared_ptr<model::Model> m_model;
    io::writer::ResultWriters     m_writer;
    mutable io::dsl::Registry     m_registry;

    std::unique_ptr<loadcase::LoadCase> m_active_loadcase;
    std::string                         m_active_loadcase_type;
    int                                 m_next_loadcase_id = 1;

public:
    Parser();
    virtual ~Parser();

    void run(const std::string&                   input_path,
             const std::string&                   output_path,
             const io::writer::WriterFileFormats& writer_formats = io::writer::WriterFileFormats());
    void document(const DocOptions& opts) const;

    model::Model& model();
    const model::Model& model() const;
    io::writer::ResultWriters& writer();
    const io::writer::ResultWriters& writer() const;
    io::dsl::Registry& registry();
    const io::dsl::Registry& registry() const;

    int next_loadcase_id();
    void set_active_loadcase(std::unique_ptr<loadcase::LoadCase> lc, std::string type);
    loadcase::LoadCase* active_loadcase();
    const loadcase::LoadCase* active_loadcase() const;
    template<class T> T* active_loadcase_as();
    template<class T> const T* active_loadcase_as() const;
    void clear_active_loadcase();
    const std::string& active_loadcase_type() const;

protected:
    virtual void configure_definition_stage(io::dsl::Registry& registry);
    virtual void configure_topology_stage  (io::dsl::Registry& registry);
    virtual void configure_field_stage     (io::dsl::Registry& registry);
    virtual void configure_data_stage      (io::dsl::Registry& registry);

private:
    void run_definition_stage(const std::string& input_path);
    void run_topology_stage  (const std::string& input_path);
    void run_field_stage     (const std::string& input_path);
    void run_data_stage      (const std::string&                   input_path,
                              const std::string&                   output_path,
                              const io::writer::WriterFileFormats& writer_formats);

    void register_topology_commands(io::dsl::Registry& registry);
    void register_analysis_commands(io::dsl::Registry& registry);
    void register_documentation_commands();
};

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
