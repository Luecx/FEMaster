#include "registry.h"
#include "command_registrar.h"
#include <stdexcept>
#include <sstream>
#include <map>
#include <string>
#include <set>
#include <iomanip>
#include <cctype>

#include "pattern.h"

namespace fem::reader2 {

Registry& Registry::instance()
{
    static Registry registry;
    return registry;
}

std::string Registry::key(const Scope& scope, const std::string& name)
{
    return scope + std::string("\x1f") + name;
}

void Registry::register_command(Scope scope, std::string name, Configure cfg) {
    CommandSpec spec;
    spec.scope = std::move(scope);
    spec.name  = std::move(name);
    if (!cfg) throw std::runtime_error("Null configure for command: " + spec.scope + ":" + spec.name);
    cfg(spec);
    const std::string compound_key = key(spec.scope, spec.name);
    _map[compound_key] = std::move(spec);
    auto& stored_spec = _map[compound_key];
    _children[stored_spec.scope].insert(stored_spec.name);
}

const CommandSpec* Registry::find(const Scope& scope, const std::string& name) const {
    auto it = _map.find(key(scope, name));
    return it == _map.end() ? nullptr : &it->second;
}

bool Registry::scope_has_children(const Scope& scope) const {
    auto it = _children.find(scope);
    return it != _children.end() && !it->second.empty();
}
bool Registry::scope_allows(const Scope& scope, const std::string& child_name) const {
    auto it = _children.find(scope);
    return it != _children.end() && it->second.count(child_name);
}
std::vector<const CommandSpec*> Registry::children_of(const Scope& scope_name) const {
    std::vector<const CommandSpec*> out;
    auto it = _children.find(scope_name);
    if (it == _children.end()) return out;
    for (const auto& child : it->second)
        if (auto* spec = find(scope_name, child)) out.push_back(spec);
    return out;
}

//--- docs helpers
static void bar(std::ostringstream& out, std::string_view s = "")
{
    out << "|  " << s << "\n";
}

static std::string render_columns(const std::vector<std::string>& names, const std::string& prefix)
{
    std::ostringstream os;
    // merge sequences like N1..N8
    auto split = [](const std::string& s, std::string& head, int& num)->bool{
        size_t pos = s.size();
        while (pos>0 && std::isdigit(static_cast<unsigned char>(s[pos-1]))) --pos;
        if (pos==s.size() || pos==0) return false;
        head = s.substr(0,pos);
        try { num = std::stoi(s.substr(pos)); } catch (...) { return false; }
        return true;
    };
    for (size_t i=0;i<names.size();){
        std::string h1; int n1=0;
        bool cr = split(names[i], h1, n1);
        size_t j = i+1;
        if (cr){
            int last = n1;
            while (j<names.size()){
                std::string hj; int nj=0;
                if (!split(names[j], hj, nj)) break;
                if (hj!=h1 || nj!=last+1) break;
                last = nj; ++j;
            }
            if (j>i+1){
                os << prefix << "[" << std::setw(2) << std::right << (i+1) << "] "
                   << std::left << std::setw(12)
                   << (h1 + std::to_string(n1) + ".." + h1 + std::to_string(n1+int(j-i-1))) << "\n";
                i = j; continue;
            }
        }
        os << prefix << "[" << std::setw(2) << std::right << (i+1) << "] "
           << std::left << std::setw(12) << names[i] << "\n";
        ++i;
    }
    return os.str();
}

//--- docs renderer
std::string Registry::document(const Scope& scope, const std::string& name) const {
    const auto* cs = find(scope, name);
    if (!cs) throw std::runtime_error("No such command: " + scope + "::" + name);

    std::ostringstream out;
    out << scope << "::" << name << "\n";
    bar(out, std::string("Title: ") + (cs->doc_title.empty() ? "(no title)" : cs->doc_title));
    if (!cs->doc_summary.empty()) {
        std::istringstream iss(cs->doc_summary);
        std::string line;
        while (std::getline(iss, line)) bar(out, "  " + line);
    }

    bar(out);
    bar(out, "Parent Scope : " + scope);
    if (cs->open_scope) bar(out, "Opens Scope  : " + *cs->open_scope);

    if (cs->key_rules) {
        bar(out);
        bar(out, "Keyword Arguments:");
        for (const auto& [key, rule] : cs->key_rules->rules()) {
            std::ostringstream head;
            head << key << (rule.required ? " (required)" : " (optional)");
            if (rule.default_value) head << ", default=" << *rule.default_value;
            bar(out, "  " + head.str());
        }
    }

    // segments
    int idx = 1;
    for (const auto& p : cs->plans) {
        std::string cond_text = p.cond().text.empty() ? "always" : p.cond().text;
        for (const auto& seg : p.segments()) {
            bar(out);
            std::ostringstream head;
            head << std::setw(2) << std::right << idx++ << ") "
                 << (seg.label().empty() ? "Segment" : seg.label())
                 << " (only possible iff " << cond_text << ")";
            bar(out, "  " + head.str());

            std::ostringstream lr;
            lr << "Groups: " << seg.range().min() << " – "
               << (seg.range().max()==kInf ? "∞" : std::to_string(seg.range().max()))
               << " | " << (seg.pattern().multiline() ? "multiline" : "single-line");
            bar(out, "    " + lr.str());

            if (!seg.pattern().doc_summary().empty())
                bar(out, "    " + seg.pattern().doc_summary());
            if (!seg.pattern().doc_notes().empty())
                bar(out, "    Notes: " + seg.pattern().doc_notes());

            bar(out, "    Columns:");
            auto names = seg.pattern().column_names();
            out << render_columns(names, "|        ");
        }
    }
    return out.str();
}

std::string Registry::document_all() const {
    std::ostringstream os;
    std::map<std::string, std::vector<std::string>> grouped;
    for (auto& [key, cs] : _map) grouped[cs.scope].push_back(cs.name);
    for (auto& [scope, cmds] : grouped) {
        os << "=== Scope " << scope << " ===\n";
        for (auto& name : cmds) os << document(scope, name) << "\n";
    }
    return os.str();
}

void Registry::print_tree(std::ostream& os) const {
    std::map<std::string, std::set<std::string>> tree;
    for (const auto& [key_str, spec] : _map) tree[spec.scope].insert(spec.name);
    std::function<void(const std::string&, int)> recurse =
        [&](const std::string& scope, int depth) {
            if (depth == 0) os << "ROOT\n";
            auto it = tree.find(scope);
            if (it == tree.end()) return;
            for (auto& child : it->second) {
                os << std::string(depth * 3, ' ') << "- " << child << "\n";
                recurse(child, depth + 1);
            }
        };
    recurse("ROOT", 0);
    os << std::endl;
}

CommandRegistrar::CommandRegistrar(const char* scope, const char* name, Registry::Configure cfg) {
    Registry::instance().register_command(scope, name, std::move(cfg));
}

} // namespace fem::reader2
