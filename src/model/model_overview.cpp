//
// Model overview pretty-printer
//

#include "model_overview.h"

#include <map>
#include <sstream>

#include "../core/logging.h"
#include "element/element.h"
#include "beam/b33.h"
#include "truss/truss.h"
#include "shell/s3.h"
#include "shell/s4.h"
#include "shell/s4_mitc.h"
#include "shell/s6.h"
#include "shell/s8.h"
#include "solid/c3d4.h"
#include "solid/c3d6.h"
#include "solid/c3d8.h"
#include "solid/c3d10.h"
#include "solid/c3d13.h"
#include "solid/c3d15.h"
#include "solid/c3d20.h"
#include "solid/c3d20r.h"
#include "../section/section_beam.h"
#include "../section/section_solid.h"
#include "../section/section_shell.h"

namespace fem { namespace model {

static std::string element_type_of(ElementInterface* e) {
    if (!e) return {};
    return e->type_name();
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

    // Top counts
    const auto n_nodes = d.node_sets.has_all() && d.node_sets.all() ? (int)d.node_sets.all()->size() : 0;
    const auto n_elems = d.elem_sets.has_all() && d.elem_sets.all() ? (int)d.elem_sets.all()->size() : 0;

    logging::info(true, "Nodes (NALL: ", n_nodes, " / max: ", (int)d.max_nodes, ")");
    logging::info(true, "Elements (EALL: ", n_elems, " / max: ", (int)d.max_elems, ")");

    // Element type distribution (best-effort via dynamic casts)
    std::map<std::string, int> by_type;
    int n_known = 0;
    for (const auto& ep : d.elements) {
        if (!ep) continue;
        auto name = element_type_of(ep.get());
        if (name.empty()) name = "UNKNOWN";
        by_type[name]++;
        ++n_known;
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

    // Sets
    print_sets_header_and_items("Node Sets",    d.node_sets);
    print_sets_header_and_items("Element Sets", d.elem_sets);
    print_sets_header_and_items("Surface Sets", d.surface_sets);
    print_sets_header_and_items("Line Sets",    d.line_sets);

    // Materials
    int n_materials = 0; for (auto it = d.materials.begin(); it != d.materials.end(); ++it) ++n_materials;
    logging::info(true, "Materials (", n_materials, ")");
    logging::up();
    for (const auto& kv : d.materials) {
        const auto& mat = kv.second;
        logging::info(true, mat ? mat->name : std::string("-"));
    }
    logging::down();
    logging::info(true, "");

    // Profiles
    int n_profiles = 0; for (auto it = d.profiles.begin(); it != d.profiles.end(); ++it) ++n_profiles;
    logging::info(true, "Profiles (", n_profiles, ")");
    logging::up();
    for (const auto& kv : d.profiles) {
        const auto& pr = kv.second;
        logging::info(true, pr ? pr->name : std::string("-"));
    }
    logging::down();
    logging::info(true, "");

    // Sections
    logging::info(true, "Sections (", (int)d.sections.size(), ")");
    logging::up();
    for (auto& s : d.sections) {
        logging::info(true, s ? s->str() : std::string("Section: (null)"));
    }
    logging::down();
    logging::info(true, "");

    // Couplings
    int n_cpl = (int)d.couplings.size();
    int n_kin = 0, n_str = 0;
    for (auto& c : d.couplings) {
        if (c.type == constraint::CouplingType::KINEMATIC) ++n_kin; else ++n_str;
    }
    logging::info(true, "Couplings (", n_cpl, ")");
    logging::up();
    logging::info(true, "Kinematic (", n_kin, ")");
    logging::up();
    for (auto& c : d.couplings) if (c.type == constraint::CouplingType::KINEMATIC) logging::info(true, c.str());
    logging::down();
    logging::info(true, "Structural (", n_str, ")");
    logging::up();
    for (auto& c : d.couplings) if (c.type == constraint::CouplingType::STRUCTURAL) logging::info(true, c.str());
    logging::down();
    logging::down();

    // Support collectors
    int n_supp_cols = 0; for (auto it = d.supp_cols.begin(); it != d.supp_cols.end(); ++it) ++n_supp_cols;
    logging::info(true, "Support Collectors (", n_supp_cols, ")");
    logging::up();
    for (const auto& kv : d.supp_cols) {
        const auto& name = kv.first; const auto& ptr = kv.second; if (!ptr) continue;
        logging::info(true, name, " (", (int)ptr->size(), ")");
        logging::up();
        for (const auto& sup : ptr->entries()) {
            logging::info(true, sup.str());
        }
        logging::down();
    }
    logging::down();

    // Load collectors
    int n_load_cols = 0; for (auto it = d.load_cols.begin(); it != d.load_cols.end(); ++it) ++n_load_cols;
    logging::info(true, "Load Collectors (", n_load_cols, ")");
    logging::up();
    for (const auto& kv : d.load_cols) {
        const auto& name = kv.first; const auto& ptr = kv.second; if (!ptr) continue;
        logging::info(true, name, " (", (int)ptr->size(), ")");
        logging::up();
        for (const auto& load_ptr : ptr->entries()) {
            if (load_ptr) logging::info(true, load_ptr->str());
        }
        logging::down();
    }
    logging::down();
}

} } // namespace fem::model
