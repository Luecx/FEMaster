//
// Model overview pretty-printer
//

#include "model_overview.h"

#include <map>
#include <sstream>

#include "../core/logging.h"
#include "element/element.h"
#include "../section/section_beam.h"
#include "../section/section_solid.h"
#include "../section/section_shell.h"

namespace fem { namespace model {

static std::string element_type_of(ElementInterface* e) {
    return e ? e->type_name() : std::string{};
}

template<typename SetsT>
static void print_sets_header_and_items(const char* title, const SetsT& sets) {
    int n_sets = 0;
    for (auto it = sets.begin(); it != sets.end(); ++it) ++n_sets;
    logging::info(true, std::string(title) + " (" + std::to_string(n_sets) + ")");
    logging::up();
    for (const auto& kv : sets) {
        const auto& name = kv.first;
        const auto& ptr  = kv.second;
        if (!ptr) continue;
        logging::info(true, "", name, " (", ptr->size(), ")");
    }
    logging::down();
    logging::info(true, "");
}

void print_model_overview(const Model& mdl) {
    const auto& d = *mdl._data;

    const Index n_nodes = d.positions ? d.positions->rows : Index(0);
    const Index n_elems = static_cast<Index>(d.elements.size());

    logging::info(true, "Nodes (", n_nodes, ")");
    logging::info(true, "Elements (", n_elems, ")");

    std::map<std::string, int> by_type;
    for (const auto& ep : d.elements) {
        if (!ep) continue;
        auto name = element_type_of(ep.get());
        if (name.empty()) name = "UNKNOWN";
        ++by_type[name];
    }

    logging::up();
    if (!by_type.empty()) {
        std::ostringstream os;
        os << "Element Types:";
        bool first = true;
        for (const auto& kv : by_type) {
            if (!first) os << "; ";
            os << " " << kv.first << " (" << kv.second << ")";
            first = false;
        }
        logging::info(true, os.str());
    } else {
        logging::info(true, "Element Types: not distinguished");
    }
    logging::down();
    logging::info(true, "");

    print_sets_header_and_items("Node Sets",    d.node_sets);
    print_sets_header_and_items("Element Sets", d.elem_sets);
    print_sets_header_and_items("Surface Sets", d.surface_sets);
    print_sets_header_and_items("Line Sets",    d.line_sets);

    int n_materials = 0;
    for (auto it = d.materials.begin(); it != d.materials.end(); ++it) ++n_materials;
    logging::info(true, "Materials (", n_materials, ")");
    logging::up();
    for (const auto& kv : d.materials) {
        const auto& mat = kv.second;
        logging::info(true, mat ? mat->name : std::string("-"));
    }
    logging::down();
    logging::info(true, "");

    int n_profiles = 0;
    for (auto it = d.profiles.begin(); it != d.profiles.end(); ++it) ++n_profiles;
    logging::info(true, "Profiles (", n_profiles, ")");
    logging::up();
    for (const auto& kv : d.profiles) {
        const auto& pr = kv.second;
        logging::info(true, pr ? pr->name : std::string("-"));
    }
    logging::down();
    logging::info(true, "");

    logging::info(true, "Sections (", static_cast<int>(d.sections.size()), ")");
    logging::up();
    for (auto& s : d.sections) {
        logging::info(true, s ? s->str() : std::string("Section: (null)"));
    }
    logging::down();
    logging::info(true, "");

    const int n_cpl = static_cast<int>(d.couplings.size());
    int n_kin = 0;
    int n_str = 0;
    for (auto& c : d.couplings) {
        if (c.type == constraint::CouplingType::KINEMATIC) ++n_kin;
        else ++n_str;
    }
    logging::info(true, "Couplings (", n_cpl, ")");
    logging::up();
    logging::info(true, "Kinematic (", n_kin, ")");
    logging::up();
    for (auto& c : d.couplings) {
        if (c.type == constraint::CouplingType::KINEMATIC) logging::info(true, c.str());
    }
    logging::down();
    logging::info(true, "Structural (", n_str, ")");
    logging::up();
    for (auto& c : d.couplings) {
        if (c.type == constraint::CouplingType::STRUCTURAL) logging::info(true, c.str());
    }
    logging::down();
    logging::down();

    int n_supp_cols = 0;
    for (auto it = d.supp_cols.begin(); it != d.supp_cols.end(); ++it) ++n_supp_cols;
    logging::info(true, "Support Collectors (", n_supp_cols, ")");
    logging::up();
    for (const auto& kv : d.supp_cols) {
        const auto& name = kv.first;
        const auto& ptr  = kv.second;
        if (!ptr) continue;
        logging::info(true, name, " (", static_cast<int>(ptr->size()), ")");
        logging::up();
        for (const auto& sup : ptr->entries()) {
            logging::info(true, sup.str());
        }
        logging::down();
    }
    logging::down();

    int n_load_cols = 0;
    for (auto it = d.load_cols.begin(); it != d.load_cols.end(); ++it) ++n_load_cols;
    logging::info(true, "Load Collectors (", n_load_cols, ")");
    logging::up();
    for (const auto& kv : d.load_cols) {
        const auto& name = kv.first;
        const auto& ptr  = kv.second;
        if (!ptr) continue;
        logging::info(true, name, " (", static_cast<int>(ptr->size()), ")");
        logging::up();
        for (const auto& load_ptr : ptr->entries()) {
            if (load_ptr) logging::info(true, load_ptr->str());
        }
        logging::down();
    }
    logging::down();
}

} } // namespace fem::model
