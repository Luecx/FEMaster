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
#include "../core/logging.h"   // fem::logging::info(...)

namespace fem {
namespace dsl {

/**
 * @class Registry
 * @brief Owns and exposes all command specifications of the DSL.
 */
struct Registry {
    std::unordered_map<std::string, Command> _map;

    /**
     * @brief Create-or-get a command by name and configure it via `fn`.
     *
     * If a command with `name` does not yet exist, it is value-initialized as
     * `Command{name}` and inserted. The provided callback `fn` is then invoked
     * with a reference to the stored `Command`, allowing in-place configuration
     * (admission rules, variants, segments, etc.).
     *
     * @param name Command name (e.g. "ELASTIC").
     * @param fn   Configuration callback receiving `Command&`.
     */
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

    /**
     * @brief Prints a compact reference of all (or one) commands.
     *
     * @param filter  Optional exact command name to restrict output; empty = all.
     * @param compact If true, per-token doc rows are hidden (short output).
     *
     * Output schema (keine '|' Pipes, nur Leerzeichen für Einrückungen):
     * @code
     * = FEM DSL — Kurzreferenz =
     *
     * *COMMAND
     * Description:
     *   ...
     * Admitted under:
     *   ...
     * Variants:
     *   - Variant #1 — When: ...
     *     Description:
     *       <variant-doc or '-'>
     *     Data-Layout:
     *       • Lines [min..max] (multiline/single-line):
     *         h1, h2, ..., g
     *         <token-name>: <doc> [<type>]
     * @endcode
     */
    void print_help(const std::string& filter = {}, bool compact = false) const {
        // Condition -> string
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
                    for (size_t i = 0; i < c.values.size(); ++i) { if (i) os << ", "; os << c.values[i]; }
                    os << " }"; return os.str();
                }
                case K::ParentKeyEquals: {
                    std::ostringstream os; os << parent << ".keys[\"" << c.key << "\"] in { ";
                    for (size_t i = 0; i < c.values.size(); ++i) { if (i) os << ", "; os << c.values[i]; }
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
                    for (size_t i = 0; i < c.values.size(); ++i) { if (i) os << ", "; os << c.values[i]; }
                    os << " }"; return os.str();
                }
                case K::All: {
                    std::ostringstream os;
                    for (size_t i = 0; i < c.children.size(); ++i) { if (i) os << " AND "; os << recurse(c.children[i], self, parent, recurse); }
                    return os.str();
                }
                case K::Any: {
                    std::ostringstream os;
                    for (size_t i = 0; i < c.children.size(); ++i) { if (i) os << " OR "; os << recurse(c.children[i], self, parent, recurse); }
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

        // Pattern render (header + optional token docs), using fem::logging::info
        auto render_pattern = [compact](const Pattern& pat) {
            std::vector<std::string> header_items;
            header_items.reserve(pat.required_tokens());
            std::size_t unnamed_next = 1;

            struct DocRow { std::string left; std::string doc; std::string type; bool is_range; };
            std::vector<DocRow> docs;

            for (const auto& pe : pat._elems) {
                const std::size_t cnt  = pe->count();
                const std::string base = pe->name_base();
                const std::string doc  = pe->description();
                const char* tname      = pe->type_name();

                if (!base.empty()) {
                    if (cnt > 1) {
                        for (std::size_t i = 1; i <= cnt; ++i) {
                            std::ostringstream h; h << base << i;
                            header_items.push_back(h.str());
                        }
                        std::ostringstream left; left << base << "1 - " << base << cnt;
                        docs.push_back({ left.str(), doc, tname, true });
                    } else {
                        header_items.push_back(base);
                        docs.push_back({ base, doc, tname, false });
                    }
                } else {
                    if (cnt > 0) {
                        for (std::size_t i = 0; i < cnt; ++i) {
                            std::ostringstream nm; nm << "t" << (unnamed_next++);
                            header_items.push_back(nm.str());
                            docs.push_back({ nm.str(), doc, tname, false });
                        }
                    }
                }
            }

            // Header line: "        h1, h2, ..., g"
            if (!header_items.empty()) {
                std::ostringstream joined;
                for (std::size_t i = 0; i < header_items.size(); ++i) {
                    if (i) joined << ", ";
                    joined << header_items[i];
                }
                fem::logging::info(true, "        ", joined.str());
            }

            if (compact) return;

            // Doc lines
            for (const auto& r : docs) {
                fem::logging::info(true,
                    "        ",
                    std::left, std::setw(14), (r.left + ":"),
                    " ", (r.doc.empty() ? "-" : r.doc), " [", r.type, "]"
                );
            }
        };

        // Sorted command names
        std::vector<std::string> names;
        names.reserve(_map.size());
        for (const auto& kv : _map) names.push_back(kv.first);
        std::sort(names.begin(), names.end());

        fem::logging::info(true, "= FEM DSL — Kurzreferenz =\n");

        for (const auto& cmd_name : names) {
            if (!filter.empty() && cmd_name != filter) continue;
            const Command& c = _map.at(cmd_name);

            fem::logging::info(true, "*", cmd_name);
            fem::logging::info(true, "Description:");
            fem::logging::info(true, "  ", (c.doc_.empty() ? std::string("-") : c.doc_));

            if (c.has_keyword_spec_) {
                const auto& entries = c.keyword_spec_.entries();
                if (!entries.empty()) {
                    std::vector<const KeywordSpec::Entry*> keyword_entries;
                    keyword_entries.reserve(entries.size());
                    for (const auto& kv : entries) {
                        keyword_entries.push_back(&kv.second);
                    }
                    std::sort(keyword_entries.begin(), keyword_entries.end(),
                        [](const KeywordSpec::Entry* a, const KeywordSpec::Entry* b) {
                            return a->canonical < b->canonical;
                        });

                    fem::logging::info(true, "Keywords:");
                    for (const auto* entry : keyword_entries) {
                        std::ostringstream header;
                        header << "  - " << entry->canonical;
                        if (entry->is_flag) {
                            header << " (flag)";
                        } else {
                            header << " (" << (entry->required ? "required" : "optional") << ")";
                        }
                        if (entry->has_default) {
                            header << ", default=\"" << entry->default_value << "\"";
                        }
                        if (!entry->allowed.empty()) {
                            header << ", allowed={";
                            for (std::size_t i = 0; i < entry->allowed.size(); ++i) {
                                if (i) header << ", ";
                                header << entry->allowed[i];
                            }
                            header << "}";
                        }
                        fem::logging::info(true, header.str());

                        const std::string doc_line = entry->doc.empty() ? std::string("-") : entry->doc;
                        fem::logging::info(true, "    ", doc_line);

                        if (!entry->alternatives.empty()) {
                            std::ostringstream aliases;
                            aliases << "    aliases: ";
                            for (std::size_t i = 0; i < entry->alternatives.size(); ++i) {
                                if (i) aliases << ", ";
                                aliases << entry->alternatives[i];
                            }
                            fem::logging::info(true, aliases.str());
                        }
                    }

                    fem::logging::info(true);
                }
            }

            fem::logging::info(true, "Admitted under:");
            fem::logging::info(true, "  ", cond_to_string(c.admit_, "self", "parent", cond_to_string));

            if (c.variants_.empty()) {
                fem::logging::info(true, "Variants:");
                fem::logging::info(true, "  (none)");
                fem::logging::info(true); // blank line
                continue;
            }

            fem::logging::info(true, "Variants:");
            for (std::size_t vi = 0; vi < c.variants_.size(); ++vi) {
                const Variant& v = c.variants_[vi];

                fem::logging::info(true,
                    "  - Variant #", (vi + 1), " — When: ",
                    (v.has_condition_
                        ? cond_to_string(v.condition_, "self", "parent", cond_to_string)
                        : std::string("(none)"))
                );

                fem::logging::info(true, "    Description:");
                fem::logging::info(true, "      ", (v.doc_.empty() ? std::string("-") : v.doc_));
                fem::logging::info(true, "    Data-Layout:");

                for (std::size_t si = 0; si < v._segments.size(); ++si) {
                    const auto& s   = v._segments[si];
                    const bool multi = s._pattern.is_multiline();

                    fem::logging::info(true,
                        "      • Lines ", fmt_range(s._range.min_, s._range.max_),
                        " (", (multi ? "multiline" : "single-line"), "):"
                    );

                    render_pattern(s._pattern);
                }
            }

            fem::logging::info(true); // blank line between commands
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
