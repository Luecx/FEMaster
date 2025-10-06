// parser/registry.cpp
/**
 * @file registry.cpp
 * @brief Implementation of the command registry and documentation utilities.
 */
#include "registry.h"
#include <stdexcept>
#include <sstream>
#include <map>
#include <string>
#include <set>

namespace fem::reader2 {

Registry& Registry::instance() {
    static Registry R;
    return R;
}

std::string Registry::key(const Scope& scope, const std::string& name) {
    return scope + std::string("\x1f") + name;
}

void Registry::register_command(Scope scope, std::string name, Configure cfg) {
    CommandSpec spec;
    spec.scope = std::move(scope);
    spec.name  = std::move(name);
    if (!cfg)
        throw std::runtime_error("Null configure for command: " + spec.scope + ":" + spec.name);
    cfg(spec);
    map_[key(spec.scope, spec.name)] = std::move(spec);
    children_[spec.scope].insert(map_[key(spec.scope, spec.name)].name);
}

const CommandSpec* Registry::find(const Scope& scope, const std::string& name) const {
    auto it = map_.find(key(scope, name));
    return it == map_.end() ? nullptr : &it->second;
}

bool Registry::scope_has_children(const Scope& scope) const {
    auto it = children_.find(scope);
    return it != children_.end() && !it->second.empty();
}

bool Registry::scope_allows(const Scope& scope, const std::string& child_name) const {
    auto it = children_.find(scope);
    return it != children_.end() && it->second.count(child_name);
}

std::vector<const CommandSpec*> Registry::children_of(const Scope& scope_name) const {
    std::vector<const CommandSpec*> out;
    auto it = children_.find(scope_name);
    if (it == children_.end()) return out;
    for (const auto& child : it->second)
        if (auto* spec = find(scope_name, child)) out.push_back(spec);
    return out;
}

std::string Registry::document(const Scope& scope, const std::string& name) const {
    auto* cs = find(scope, name);
    if (!cs)
        throw std::runtime_error("No such command: " + scope + ":" + name);

    std::ostringstream md;
    md << "### `" << scope << "::" << name << "`\n";
    md << "- Parent: `" << scope << "`\n";
    md << "- Has data lines: " << (cs->has_lines ? "yes" : "no") << "\n";
    if (cs->open_scope) {
        md << "- Opens scope: `" << *cs->open_scope << "`\n";
    }
    md << "\n";
    return md.str();
}

std::string Registry::document_all() const {
    std::map<std::string, std::vector<std::string>> grouped;
    for (auto& [k, cs] : map_) grouped[cs.scope].push_back(cs.name);
    std::ostringstream md;
    for (auto& [scope, cmds] : grouped) {
        md << "## Scope `" << scope << "`\n";
        for (auto& c : cmds) md << document(scope, c);
    }
    return md.str();
}

CommandRegistrar::CommandRegistrar(const char* scope, const char* name, Registry::Configure cfg) {
    Registry::instance().register_command(scope, name, std::move(cfg));
}

void Registry::print_tree(std::ostream& os) const {
    // Build reverse map: parent_scope -> child names
    std::map<std::string, std::set<std::string>> tree;
    for (const auto& [key_str, spec] : map_) {
        tree[spec.scope].insert(spec.name);
    }

    // Recursive printer
    std::function<void(const std::string&, int)> recurse =
        [&](const std::string& scope, int depth) {
            if (depth == 0) os << "ROOT\n";
            auto it = tree.find(scope);
            if (it == tree.end()) return;

            for (auto& child : it->second) {
                os << std::string(depth * 3, ' ') << "├── " << child << "\n";
                recurse(child, depth + 1);
            }
    };

    recurse("ROOT", 0);
    os << std::endl;
}

} // namespace fem::reader2
