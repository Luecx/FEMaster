#pragma once

#include <memory>
#include <string>

#include "../dsl/registry.h"
#include "../reader/writer.h"

namespace fem {
namespace loadcase { struct LoadCase; }
namespace dsl      { class File; struct Line; }
namespace model    { struct Model; }

namespace input_decks {

// ---- Doc mode options ----
struct DocOptions {
  enum class Action { List, Show, Tokens, Variants, Search, WhereToken, All };
  enum class Format { Text, Markdown, Json };
  enum class Verbosity { Index, Compact, Full };

  Action action = Action::List;
  std::string cmd;         // for Show/Tokens/Variants
  std::string query;       // for Search / WhereToken
  Format format = Format::Text;
  Verbosity verbosity = Verbosity::Full;
  int wrap_width = 100;    // 0=no wrap
  bool regex = false;
  bool no_wrap = false;
};

class Parser {
public:
    // Construct with a tiny placeholder model and register all commands immediately.
    Parser();
    ~Parser();

    // High-level modes
    void run(const std::string& input_path, const std::string& output_path);
    void document(const DocOptions& opts) const;

    // Accessors
    model::Model& model();
    const model::Model& model() const;
    reader::Writer& writer();
    const reader::Writer& writer() const;
    dsl::Registry& registry();
    const dsl::Registry& registry() const;

    // Loadcase helpers used by command registrations
    int next_loadcase_id();
    void set_active_loadcase(std::unique_ptr<loadcase::LoadCase> lc, std::string type);
    loadcase::LoadCase* active_loadcase();
    const loadcase::LoadCase* active_loadcase() const;
    template<class T> T* active_loadcase_as();
    template<class T> const T* active_loadcase_as() const;
    void clear_active_loadcase();
    const std::string& active_loadcase_type() const;

private:
    struct AnalyseData {
        int highest_node_id    = -1;
        int highest_element_id = -1;
        int highest_surface_id = -1;
    };

    // Pipeline steps
    AnalyseData preprocess_ids(const std::string& input_path);
    static std::size_t required_nodes_for_type(const std::string& type);

    void analyse_nodes(dsl::File& file, dsl::Line& current, AnalyseData& A);
    void analyse_elements(dsl::File& file, dsl::Line& current, AnalyseData& A);
    void analyse_surfaces(dsl::File& file, dsl::Line& current, AnalyseData& A);

    // (Re)build model and (re)register commands
    void build_model_from(const AnalyseData& A);
    void register_commands(); // uses current m_model & resets m_registry if needed

private:
    // Runtime state
    std::shared_ptr<model::Model> m_model;
    reader::Writer                m_writer;   // opened in run()
    mutable dsl::Registry         m_registry; // re-created when model changes

    std::unique_ptr<loadcase::LoadCase> m_active_loadcase;
    std::string                         m_active_loadcase_type;
    int                                  m_next_loadcase_id = 1;

    bool m_commands_registered = false;
};

// Templates
template<class T>
inline T* Parser::active_loadcase_as() {
    return dynamic_cast<T*>(active_loadcase());
}
template<class T>
inline const T* Parser::active_loadcase_as() const {
    return dynamic_cast<const T*>(active_loadcase());
}

} // namespace input_decks
} // namespace fem
