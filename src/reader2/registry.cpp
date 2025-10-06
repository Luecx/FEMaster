#include "registry.h"
#include "command_registrar.h"
#include <stdexcept>
#include <sstream>
#include <map>
#include <string>
#include <set>
#include <iomanip>
#include <cctype>
#include <algorithm>
#include <vector>

#include "pattern.h"

namespace fem::reader2 {

Registry& Registry::instance() {
    static Registry registry;
    return registry;
}

std::string Registry::key(const Scope& scope, const std::string& name) {
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

//------------------------------------------------------------------------------
// docs helpers
//------------------------------------------------------------------------------

static void bar(std::ostringstream& out, std::string_view s = "") {
    out << "|  " << s << "\n";
}

static void bar_indent(std::ostringstream& out, int spaces, std::string_view s = "") {
    out << "|  " << std::string(spaces, ' ') << s << "\n";
}

static std::string join_sorted(const std::unordered_set<std::string>& vals, const char* sep) {
    if (vals.empty()) return {};
    std::vector<std::string> v(vals.begin(), vals.end());
    std::sort(v.begin(), v.end());
    std::ostringstream os;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) os << sep;
        os << v[i];
    }
    return os.str();
}

static const char* kind_name(PatternElement::ValueKind k) {
    switch (k) {
        case PatternElement::ValueKind::Integer:  return "Integer";
        case PatternElement::ValueKind::Floating: return "Floating";
        case PatternElement::ValueKind::Text:     return "Text";
        default:                                  return "?";
    }
}

// Render a compact, aligned table of column specs for a pattern.
// Format per row:
//   [idx]  Name(s)  Kind  Count  : Description
static std::string render_columns_table(const Pattern& pat, const std::string& prefix) {
    struct Row {
        size_t      idx = 0;          // 1-based start column index
        std::string name;             // merged name(s), e.g., N1..N8
        std::string kind;             // Integer | Floating | Text
        std::string count;            // 1 | N | min..max (∞ allowed)
        std::string desc;             // description text
    };

    std::vector<Row> rows;
    rows.reserve(pat.elements().size());

    // Determine starting index per element using MAX counts for prior elements
    size_t idx = 1;
    for (const auto& el : pat.elements()) {
        Row r;
        r.idx  = idx;

        const size_t minc = el.min_count();
        const size_t maxc = el.max_count();

        // name(s)
        if (maxc <= 1) {
            r.name = el.base_name();
        } else {
            // base1..baseN, allow ∞
            std::ostringstream nm;
            nm << el.base_name() << "1.." << el.base_name()
               << (maxc == kInf ? "Infty" : std::to_string(maxc));
            r.name = nm.str();
        }

        // kind
        r.kind = kind_name(el.kind());

        // count
        if (el.type() == PatternElement::Type::Variable) {
            std::ostringstream cs;
            cs << minc << ".." << (maxc == kInf ? "Infty" : std::to_string(maxc));
            r.count = cs.str();
        } else if (el.type() == PatternElement::Type::Fixed) {
            r.count = std::to_string(maxc);
        } else { // Single
            r.count = "1";
        }

        // description
        r.desc = el.description();

        rows.push_back(std::move(r));

        // advance index by MAX extent
        idx += maxc;
    }

    // Compute widths for alignment
    size_t name_w  = std::string("Name(s)").size();
    size_t kind_w  = std::string("Kind").size();
    size_t count_w = std::string("Count").size();
    for (const auto& r : rows) {
        name_w  = std::max(name_w,  r.name.size());
        kind_w  = std::max(kind_w,  r.kind.size());
        count_w = std::max(count_w, r.count.size());
    }

    std::ostringstream os;

    // Header row
    os << prefix << "[idx] "
       << std::left << std::setw(static_cast<int>(name_w))  << "Name(s)" << "  "
       << std::left << std::setw(static_cast<int>(kind_w))  << "Kind"    << "  "
       << std::left << std::setw(static_cast<int>(count_w)) << "Count"   << "  : "
       << "Description" << "\n";

    // Data rows
    for (const auto& r : rows) {
        os << prefix
           << "[" << std::setw(2) << std::right << r.idx << "]  "
           << std::left << std::setw(static_cast<int>(name_w))  << r.name  << "  "
           << std::left << std::setw(static_cast<int>(kind_w))  << r.kind  << "  "
           << std::left << std::setw(static_cast<int>(count_w)) << r.count << "  : "
           << r.desc << "\n";
    }

    return os.str();
}

//------------------------------------------------------------------------------
// docs renderer
//------------------------------------------------------------------------------

std::string Registry::document(const Scope& scope, const std::string& name) const {
    const auto* cs = find(scope, name);
    if (!cs) throw std::runtime_error("No such command: " + scope + "::" + name);

    std::ostringstream out;

    // Header: *NAME or *SCOPE -> *NAME
    {
        std::string header;
        if (scope == "ROOT") {
            header = "*" + cs->name;
        } else {
            header = "*" + scope + " -> *" + cs->name;
        }
        out << header << "\n";
    }

    // Description
    {
        bar(out, "Description:");
        if (!cs->doc_summary.empty()) {
            std::istringstream iss(cs->doc_summary);
            std::string line;
            while (std::getline(iss, line)) {
                bar_indent(out, 2, line);
            }
        }
    }

    // Keyword arguments
    if (cs->key_rules) {
        bar(out);
        bar(out, "Keyword arguments:");

        // compute left label width for aligned colon
        size_t left_width = 0;
        struct Row { std::string key, left, desc, values; bool has_values; };
        std::vector<Row> rows;

        for (const auto& [key, rule] : cs->key_rules->rules()) {
            std::ostringstream left;
            left << key << (rule.required ? " (required)" : " (optional)");
            if (!rule.required && rule.default_value) {
                left << " [default=" << *rule.default_value << "]";
            }
            std::string left_s = left.str();
            left_width = std::max(left_width, left_s.size());

            Row r;
            r.key = key;
            r.left = std::move(left_s);
            r.desc = rule.description;
            if (!rule.allowed_values.empty()) {
                r.values = join_sorted(rule.allowed_values, ", ");
                r.has_values = true;
            } else {
                r.values = "ANY";
                r.has_values = true;
            }
            rows.push_back(std::move(r));
        }

        // render rows
        for (const auto& r : rows) {
            std::ostringstream line;
            line << std::left << std::setw(static_cast<int>(left_width)) << r.left << " : ";
            if (!r.desc.empty()) line << r.desc;
            if (r.has_values) {
                if (!r.desc.empty()) line << " ";
                line << "[possible values: " << r.values << "]";
            }
            bar_indent(out, 2, line.str());
        }
    }

    // Variants
    {
        bar(out);
        bar(out, "Variants:");
        int idx = 1;
        for (const auto& p : cs->plans) {
            std::string cond_text = p.cond().text.empty() ? "always" : p.cond().text;
            for (const auto& seg : p.segments()) {
                bar_indent(out, 2); // empty spacer line between variants
                {
                    std::ostringstream head;
                    head << "(" << idx++ << ") "
                         << (seg.label().empty() ? "Segment" : seg.label());
                    bar_indent(out, 2, head.str());
                }
                bar_indent(out, 6, std::string("Condition    : ") + cond_text);

                std::ostringstream lr;
                lr << "Groups       : " << seg.range().min() << " – "
                   << (seg.range().max()==kInf ? "Infty" : std::to_string(seg.range().max()))
                   << " | " << (seg.pattern().multiline() ? "multiline" : "single-line");
                bar_indent(out, 6, lr.str());

                if (!seg.pattern().doc_summary().empty())
                    bar_indent(out, 6, std::string("Summary      : ") + seg.pattern().doc_summary());
                if (!seg.pattern().doc_notes().empty())
                    bar_indent(out, 6, std::string("Notes        : ") + seg.pattern().doc_notes());

                bar_indent(out, 6, "Columns:");
                out << render_columns_table(seg.pattern(), "|            ");
            }
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

} // namespace fem::reader2
