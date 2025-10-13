#pragma once

#include <memory>
#include <string>

#include "../dsl/registry.h"
#include "../reader/writer.h"

namespace fem {
namespace loadcase {
struct LoadCase;
} // namespace loadcase

namespace dsl {
class File;
struct Line;
}

namespace model { struct Model; }

namespace input_decks {

class Parser {
public:
    Parser(std::string input_path, std::string output_path);
    ~Parser();

    void preprocess();
    void parse();

    model::Model& model();
    const model::Model& model() const;
    reader::Writer& writer();
    const reader::Writer& writer() const;
    dsl::Registry& registry();

    void register_all_commands();
    void prepare_for_documentation();

    // Loadcase helpers used by command registrations
    int next_loadcase_id();
    void set_active_loadcase(std::unique_ptr<loadcase::LoadCase> lc, std::string type);
    loadcase::LoadCase* active_loadcase();
    const loadcase::LoadCase* active_loadcase() const;
    template<class T>
    T* active_loadcase_as();
    template<class T>
    const T* active_loadcase_as() const;
    void clear_active_loadcase();
    const std::string& active_loadcase_type() const;

private:
    struct AnalyseData {
        int highest_node_id = -1;
        int highest_element_id = -1;
        int highest_surface_id = -1;
    };

    void analyse();
    void analyse_nodes(dsl::File& file, dsl::Line& current);
    void analyse_elements(dsl::File& file, dsl::Line& current);
    void analyse_surfaces(dsl::File& file, dsl::Line& current);

    void register_commands();

    std::string m_input_path;
    std::string m_output_path;

    AnalyseData m_analyse;

    std::shared_ptr<model::Model> m_model;
    reader::Writer m_writer;
    dsl::Registry m_registry;

    std::unique_ptr<loadcase::LoadCase> m_active_loadcase;
    std::string m_active_loadcase_type;
    int m_next_loadcase_id = 1;

    bool m_commands_registered = false;
};

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
