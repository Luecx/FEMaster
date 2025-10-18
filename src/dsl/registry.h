#pragma once
#include <functional>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <regex>

#include "command.h"
#include "condition.h"
#include "pattern.h"
#include "pattern_element.h"
#include "segment.h"

namespace fem {
namespace dsl {

/**
 * @class Registry
 * @brief Owns and exposes all command specifications of the DSL.
 */
struct Registry {
    std::unordered_map<std::string, Command> _map;

    void command(const std::string& name, const std::function<void(Command&)>& fn) {
        auto it = _map.find(name);
        if (it == _map.end())
            it = _map.emplace(name, Command{name}).first;
        fn(it->second);
    }

    Command* find(const std::string& name) {
        auto it = _map.find(name);
        return it == _map.end() ? nullptr : &it->second;
    }
    const Command* find(const std::string& name) const {
        auto it = _map.find(name);
        return it == _map.end() ? nullptr : &it->second;
    }

    // ---------------------------------------------------------------------
    // Documentation helpers
    // ---------------------------------------------------------------------

    // Index (just names)
    void print_index() const {
        std::vector<std::string> names;
        names.reserve(_map.size());
        for (const auto& kv : _map) names.push_back(kv.first);
        std::sort(names.begin(), names.end());
        std::cout << "= FEM DSL — Command Index =\n\n";
        for (const auto& n : names) std::cout << n << "\n";
    }

    // Compact/full help for one or all commands (filter empty => all).
    // 'compact' hides per-token doc rows to keep output short.
    void print_help(const std::string& filter = {}, bool compact = false) const {
        auto cond_to_string = [](const Condition& c,
                                 const std::string& self,
                                 const std::string& parent,
                                 const auto& recurse) -> std::string
        {
            using K = Condition::Kind;
            switch (c.kind) {
                case K::Always: return "(none)";

                case K::ParentCommandIs: {
                    std::ostringstream os; os << parent << ".command in { ";
                    for (size_t i=0;i<c.values.size();++i) { if (i) os << ", "; os << c.values[i]; }
                    os << " }"; return os.str();
                }
                case K::ParentKeyEquals: {
                    std::ostringstream os; os << parent << ".keys[\"" << c.key << "\"] in { ";
                    for (size_t i=0;i<c.values.size();++i) { if (i) os << ", "; os << c.values[i]; }
                    os << " }"; return os.str();
                }
                case K::ParentHasKey: {
                    std::ostringstream os; os << parent << ".has_key(\"" << c.key << "\")"; return os.str();
                }
                case K::KeyPresent: {
                    std::ostringstream os; os << self << ".has_key(\"" << c.key << "\")"; return os.str();
                }
                case K::KeyEquals: {
                    std::ostringstream os; os << self << ".keys[\"" << c.key << "\"] in { ";
                    for (size_t i=0;i<c.values.size();++i) { if (i) os << ", "; os << c.values[i]; }
                    os << " }"; return os.str();
                }
                case K::All: {
                    std::ostringstream os;
                    for (size_t i=0;i<c.children.size();++i) {
                        if (i) os << " AND ";
                        os << recurse(c.children[i], self, parent, recurse);
                    }
                    return os.str();
                }
                case K::Any: {
                    std::ostringstream os;
                    for (size_t i=0;i<c.children.size();++i) {
                        if (i) os << " OR ";
                        os << recurse(c.children[i], self, parent, recurse);
                    }
                    return os.str();
                }
                case K::Not: {
                    std::ostringstream os;
                    os << "NOT(" << (c.children.empty() ? "(none)" : recurse(c.children.front(), self, parent, recurse)) << ")";
                    return os.str();
                }
            }
            return "";
        };

        auto fmt_range = [](std::size_t min_v, std::size_t max_v) -> std::string {
            std::ostringstream os;
            os << "[" << min_v << "..";
            if (max_v == std::numeric_limits<std::size_t>::max()) os << "∞";
            else os << max_v;
            os << "]";
            return os.str();
        };

        auto render_pattern = [compact](const Pattern& pat, std::ostream& os) {
            // Build header items (explicit enumeration for named multi-count elements).
            std::vector<std::string> header_items;
            header_items.reserve(pat.required_tokens());

            // For unnamed elements we still generate t1, t2, ... with a running counter.
            std::size_t unnamed_next = 1;

            struct DocRow { std::string left; std::string doc; std::string type; bool is_range; };
            std::vector<DocRow> docs;

            for (const auto& pe : pat._elems) {
                const std::size_t cnt = pe->count();
                const std::string base = pe->name_base();
                const std::string doc  = pe->description();
                const char* tname      = pe->type_name();

                if (!base.empty()) {
                    if (cnt > 1) {
                        for (std::size_t i=1; i<=cnt; ++i) {
                            std::ostringstream h; h << base << i;
                            header_items.push_back(h.str());
                        }
                        std::ostringstream left; left << base << "1 - " << base << cnt;
                        docs.push_back({ left.str(), doc, tname, true });
                    } else {
                        header_items.push_back(base);
                        std::ostringstream left; left << base;
                        docs.push_back({ left.str(), doc, tname, false });
                    }
                } else {
                    if (cnt > 0) {
                        for (std::size_t i=0; i<cnt; ++i) {
                            std::ostringstream nm; nm << "t" << (unnamed_next++);
                            header_items.push_back(nm.str());
                            std::ostringstream left; left << nm.str();
                            docs.push_back({ left.str(), doc, tname, false });
                        }
                    }
                }
            }

            // Header line (explicit names)
            os << "|         ";
            for (std::size_t i=0; i<header_items.size(); ++i) {
                if (i) os << ", ";
                os << header_items[i];
            }
            os << "\n";

            if (compact) return;

            // Doc lines
            for (const auto& r : docs) {
                os << "|         " << std::left << std::setw(12) << (r.left + ":")
                   << " " << (r.doc.empty() ? "-" : r.doc) << " [" << r.type << "]\n";
            }
        };

        // Sorted output
        std::vector<std::string> names;
        names.reserve(_map.size());
        for (const auto& kv : _map) names.push_back(kv.first);
        std::sort(names.begin(), names.end());

        std::cout << "= FEM DSL — Kurzreferenz =\n\n";

        for (const auto& cmd_name : names) {
            if (!filter.empty() && cmd_name != filter) continue;
            const Command& c = _map.at(cmd_name);

            std::cout << "*" << cmd_name << "\n";
            std::cout << "| Description:\n";
            std::cout << "|   " << (c.doc_.empty() ? "-" : c.doc_) << "\n";
            std::cout << "| Admitted under:\n";
            std::cout << "|   " << cond_to_string(c.admit_, "self", "parent", cond_to_string) << "\n";

            if (c.variants_.empty()) {
                std::cout << "| Variants:\n";
                std::cout << "|   (none)\n|\n";
                continue;
            }

            std::cout << "| Variants:\n";
            for (std::size_t vi=0; vi<c.variants_.size(); ++vi) {
                const Variant& v = c.variants_[vi];

                std::cout << "|   - Variant #" << (vi+1) << " — When: "
                          << (v.has_condition_ ? cond_to_string(v.condition_, "self", "parent", cond_to_string) : "(none)") << "\n";

                std::cout << "|     Description:\n";
                std::cout << "|       " << (v.doc_.empty() ? "-" : v.doc_) << "\n";

                std::cout << "|     Data-Layout:\n";
                for (std::size_t si=0; si<v._segments.size(); ++si) {
                    const auto& s = v._segments[si];
                    const bool multi = s._pattern.is_multiline();
                    std::cout << "|       • Lines " << fmt_range(s._range.min_, s._range.max_)
                              << " (" << (multi ? "multiline" : "single-line") << "):\n";
                    render_pattern(s._pattern, std::cout);
                }
            }
            std::cout << "|\n";
        }
    }

    void print_variants(const std::string& cmd) const {
        const Command* c = find(cmd);
        if (!c) { std::cout << "Command not found: " << cmd << "\n"; return; }

        auto cond_to_string = [](const Condition& c,
                                 const std::string& self,
                                 const std::string& parent,
                                 const auto& recurse) -> std::string
        {
            using K = Condition::Kind;
            switch (c.kind) {
                case K::Always: return "(none)";
                case K::ParentCommandIs: {
                    std::ostringstream os; os << parent << ".command in { ";
                    for (size_t i=0;i<c.values.size();++i) { if (i) os << ", "; os << c.values[i]; }
                    os << " }"; return os.str();
                }
                case K::ParentKeyEquals: {
                    std::ostringstream os; os << parent << ".keys[\"" << c.key << "\"] in { ";
                    for (size_t i=0;i<c.values.size();++i) { if (i) os << ", "; os << c.values[i]; }
                    os << " }"; return os.str();
                }
                case K::ParentHasKey: {
                    std::ostringstream os; os << parent << ".has_key(\"" << c.key << "\")"; return os.str();
                }
                case K::KeyPresent: {
                    std::ostringstream os; os << self << ".has_key(\"" << c.key << "\")"; return os.str();
                }
                case K::KeyEquals: {
                    std::ostringstream os; os << self << ".keys[\"" << c.key << "\"] in { ";
                    for (size_t i=0;i<c.values.size();++i) { if (i) os << ", "; os << c.values[i]; }
                    os << " }"; return os.str();
                }
                case K::All: {
                    std::ostringstream os;
                    for (size_t i=0;i<c.children.size();++i) {
                        if (i) os << " AND ";
                        os << recurse(c.children[i], self, parent, recurse);
                    }
                    return os.str();
                }
                case K::Any: {
                    std::ostringstream os;
                    for (size_t i=0;i<c.children.size();++i) {
                        if (i) os << " OR ";
                        os << recurse(c.children[i], self, parent, recurse);
                    }
                    return os.str();
                }
                case K::Not: {
                    std::ostringstream os;
                    os << "NOT(" << (c.children.empty() ? "(none)" : recurse(c.children.front(), self, parent, recurse)) << ")";
                    return os.str();
                }
            }
            return "";
        };

        std::cout << "= Variants — " << cmd << " =\n\n";
        if (c->variants_.empty()) { std::cout << "(none)\n"; return; }

        for (std::size_t vi=0; vi<c->variants_.size(); ++vi) {
            const Variant& v = c->variants_[vi];
            std::cout << "- Variant #" << (vi+1) << " — When: "
                      << (v.has_condition_ ? cond_to_string(v.condition_, "self", "parent", cond_to_string) : "(none)") << "\n";
            std::cout << "  Doc: " << (v.doc_.empty() ? "-" : v.doc_) << "\n";
        }
    }

    void print_tokens(const std::string& cmd) const {
        const Command* c = find(cmd);
        if (!c) { std::cout << "Command not found: " << cmd << "\n"; return; }

        auto render_pattern_compact = [](const Pattern& pat, std::ostream& os) {
            std::vector<std::string> header_items;
            header_items.reserve(pat.required_tokens());
            std::size_t unnamed_next = 1;

            struct DocRow { std::string left; std::string doc; std::string type; bool is_range; };
            std::vector<DocRow> docs;

            for (const auto& pe : pat._elems) {
                const std::size_t cnt = pe->count();
                const std::string base = pe->name_base();
                const std::string doc  = pe->description();
                const char* tname      = pe->type_name();

                if (!base.empty()) {
                    if (cnt > 1) {
                        for (std::size_t i=1; i<=cnt; ++i) {
                            std::ostringstream h; h << base << i;
                            header_items.push_back(h.str());
                        }
                        std::ostringstream left; left << base << "1 - " << base << cnt;
                        docs.push_back({ left.str(), doc, tname, true });
                    } else {
                        header_items.push_back(base);
                        std::ostringstream left; left << base;
                        docs.push_back({ left.str(), doc, tname, false });
                    }
                } else {
                    if (cnt > 0) {
                        for (std::size_t i=0; i<cnt; ++i) {
                            std::ostringstream nm; nm << "t" << (unnamed_next++);
                            header_items.push_back(nm.str());
                            std::ostringstream left; left << nm.str();
                            docs.push_back({ left.str(), doc, tname, false });
                        }
                    }
                }
            }

            os << "Header: ";
            for (std::size_t i=0; i<header_items.size(); ++i) {
                if (i) os << ", ";
                os << header_items[i];
            }
            os << "\n";

            for (const auto& r : docs) {
                os << "  - " << std::left << std::setw(12) << (r.left + ":")
                   << " " << (r.doc.empty() ? "-" : r.doc) << " [" << r.type << "]\n";
            }
        };

        std::cout << "= Tokens — " << cmd << " =\n\n";
        if (c->variants_.empty()) { std::cout << "(none)\n"; return; }

        for (std::size_t vi=0; vi<c->variants_.size(); ++vi) {
            const Variant& v = c->variants_[vi];
            std::cout << "Variant #" << (vi+1) << ":\n";
            for (std::size_t si=0; si<v._segments.size(); ++si) {
                const auto& s = v._segments[si];
                render_pattern_compact(s._pattern, std::cout);
            }
            std::cout << "\n";
        }
    }

    void print_search(const std::string& query, bool regex=false) const {
        std::cout << "= Search: " << query << " =\n\n";
        std::vector<std::string> names;
        names.reserve(_map.size());
        for (const auto& kv : _map) names.push_back(kv.first);
        std::sort(names.begin(), names.end());

        std::regex rx;
        if (regex) rx = std::regex(query, std::regex::icase);

        int hits = 0;
        for (const auto& name : names) {
            const Command& c = _map.at(name);
            const auto in = [&](const std::string& s)->bool{
                if (regex) {
                    return std::regex_search(s, rx);
                } else {
                    auto hay = s; auto needle = query;
                    std::transform(hay.begin(), hay.end(), hay.begin(), ::tolower);
                    std::transform(needle.begin(), needle.end(), needle.begin(), ::tolower);
                    return hay.find(needle) != std::string::npos;
                }
            };
            bool match = in(name) || in(c.doc_);
            if (match) {
                std::cout << name << "  —  " << (c.doc_.empty() ? "-" : c.doc_) << "\n";
                ++hits;
            }
        }
        if (hits == 0) std::cout << "(no matches)\n";
    }

    void print_where_token(const std::string& token_substr) const {
        std::cout << "= Where Token: " << token_substr << " =\n\n";
        std::vector<std::string> names;
        names.reserve(_map.size());
        for (const auto& kv : _map) names.push_back(kv.first);
        std::sort(names.begin(), names.end());

        auto contains_ci = [](const std::string& s, const std::string& sub)->bool{
            std::string a = s, b = sub;
            std::transform(a.begin(), a.end(), a.begin(), ::tolower);
            std::transform(b.begin(), b.end(), b.begin(), ::tolower);
            return a.find(b) != std::string::npos;
        };

        int hits = 0;
        for (const auto& nm : names) {
            const Command& c = _map.at(nm);
            bool found = false;

            for (const auto& v : c.variants_) {
                for (const auto& seg : v._segments) {
                    for (const auto& pe : seg._pattern._elems) {
                        if (contains_ci(pe->name_base(), token_substr) || contains_ci(pe->description(), token_substr)) {
                            found = true; break;
                        }
                    }
                    if (found) break;
                }
                if (found) break;
            }

            if (found) { std::cout << nm << "\n"; ++hits; }
        }
        if (hits == 0) std::cout << "(no matches)\n";
    }
};

} // namespace dsl
} // namespace fem
