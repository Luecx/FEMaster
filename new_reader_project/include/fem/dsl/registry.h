// ===== FILE: ./include/fem/dsl/registry.h =====
/**
 * @file
 * @brief Command registry for the DSL: declares and looks up command specifications.
 *
 * The `Registry` owns all `Command` specifications keyed by their command name
 * (e.g. `"ELASTIC"`, `"NUMEIGENVALUES"`). It provides:
 *
 *  - `command(name, fn)`: create-or-get and configure a `Command` in-place via a callback.
 *  - `find(name)`: retrieve a pointer to a registered `Command` (const and non-const).
 *  - `print_help([filter]) const`: header-only console documentation using `std::cout`.
 *
 * Notes:
 *  - Registration is idempotent per name; calling `command()` multiple times for the same
 *    name will reuse and further configure the existing `Command`.
 *  - The help printer renders compact output per Command/Variant and uses Pattern elements
 *    (`name(base)`, `desc(...)`, `count()`, `type_name()`) to label tokens.
 */

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
    /**
     * @brief Storage mapping command names to their specifications.
     */
    std::unordered_map<std::string, Command> _map;

    /**
     * @brief Create-or-get a command by name and configure it via `fn`.
     *
     * If a command with `name` does not yet exist, it is value-initialized as
     * `Command{name}` and inserted. The provided callback `fn` is then invoked
     * with a reference to the stored `Command`, allowing in-place configuration
     * (admission rules, variants, segments, etc.).
     *
     * @param name Command name (e.g. `"ELASTIC"`).
     * @param fn   Configuration callback receiving `Command&`.
     */
    void command(const std::string& name, const std::function<void(Command&)>& fn) {
        auto it = _map.find(name);
        if (it == _map.end())
            it = _map.emplace(name, Command{name}).first;
        fn(it->second);
    }

    /**
     * @brief Finds a command by name (mutable access).
     *
     * @param name Command name to look up.
     * @return Pointer to the stored `Command`, or `nullptr` if not found.
     */
    Command* find(const std::string& name) {
        auto it = _map.find(name);
        return it == _map.end() ? nullptr : &it->second;
    }

    /**
     * @brief Finds a command by name (const access).
     *
     * @param name Command name to look up.
     * @return Pointer to the stored `Command` (const), or `nullptr` if not found.
     */
    const Command* find(const std::string& name) const {
        auto it = _map.find(name);
        return it == _map.end() ? nullptr : &it->second;
    }

    // ---------------------------------------------------------------------
    // Header-only console help output (std::cout)
    // ---------------------------------------------------------------------

    /**
     * @brief Prints a compact reference of all (or one) commands to `std::cout`.
     *
     * @param filter Optional exact command name to restrict output; empty = all.
     *
     * Output schema:
     * @code
     * = FEM DSL — Kurzreferenz =
     *
     * *COMMAND
     * | Description:
     * |   ...
     * | Admitted under:
     * |   ...
     * | Variants:
     * |   - Variant #1 — When: ...
     * |     Description:
     * |       <variant-doc or '-'>
     * |     Data-Layout:
     * |       • Lines [min..max] (multiline/single-line):
     * |         header (explicit names, e.g., h1, h2, ..., h10, g)
     * |         per-element docs:
     * |           - named multi-count element as a single line (e.g. h1 - h10: ...)
     * |           - others as one line per token (or per element if single)
     * |
     * @endcode
     */
    void print_help(const std::string& filter = {}) const {
        /**
         * @brief Renders a `Condition` tree into a compact, human-readable string.
         *
         * @param c        Condition to render.
         * @param self     Symbolic name representing "self" in the docs (usually "self").
         * @param parent   Symbolic name representing "parent" in the docs (usually "parent").
         * @param recurse  Reference to this lambda for recursion.
         * @return A compact string representation of the condition.
         */
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

        /**
         * @brief Formats a line range as `[min..max]`, printing `∞` for an unbounded max.
         *
         * @param min_v Inclusive minimum line count.
         * @param max_v Inclusive maximum line count.
         * @return A formatted string for the range.
         */
        auto fmt_range = [](std::size_t min_v, std::size_t max_v) -> std::string {
            std::ostringstream os;
            os << "[" << min_v << "..";
            if (max_v == std::numeric_limits<std::size_t>::max()) os << "∞";
            else os << max_v;
            os << "]";
            return os.str();
        };

        /**
         * @brief Renders a `Pattern` into compact documentation.
         *
         * Rules:
         *  - The header line (CSV) **expands** any named multi-count element `base[N]`
         *    to explicit names: `base1, base2, ..., baseN`.
         *  - For docs:
         *      • A **named** multi-count element prints **one** line: `base1 - baseN:  desc [type]`.
         *      • A **named** single-count element prints one line: `base:  desc [type]`.
         *      • An **unnamed** element falls back to synthetic names `tK` (expanded).
         *
         * @param pat Pattern to render.
         * @param os  Output stream to write to (usually std::cout).
         */
        auto render_pattern = [](const Pattern& pat, std::ostream& os) {
            // Build header items (explicit enumeration for named multi-count elements).
            std::vector<std::string> header_items;
            header_items.reserve(pat.required_tokens());

            // For unnamed elements we still generate t1, t2, ... with a running counter.
            std::size_t unnamed_next = 1;

            // Also collect doc rows:
            struct DocRow { std::string left; std::string doc; std::string type; bool is_range; };
            std::vector<DocRow> docs;

            for (const auto& pe : pat._elems) {
                const std::size_t cnt = pe->count();
                const std::string base = pe->name_base();
                const std::string doc  = pe->description();
                const char* tname      = pe->type_name();

                if (!base.empty()) {
                    if (cnt > 1) {
                        // Header explicit: base1, base2, ..., baseN
                        for (std::size_t i=1; i<=cnt; ++i) {
                            std::ostringstream h; h << base << i;
                            header_items.push_back(h.str());
                        }
                        // Doc compact: base1 - baseN: doc [type]
                        {
                            std::ostringstream left; left << base << "1 - " << base << cnt;
                            docs.push_back({ left.str(), doc, tname, true });
                        }
                    } else {
                        // Single named element
                        header_items.push_back(base);
                        std::ostringstream left; left << base;
                        docs.push_back({ left.str(), doc, tname, false });
                    }
                } else {
                    // Unnamed element: enumerate synthetic names tK
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

            // Doc lines
            for (const auto& r : docs) {
                os << "|         " << std::left << std::setw(12) << (r.left + ":")
                   << " " << (r.doc.empty() ? "-" : r.doc) << " [" << r.type << "]\n";
            }
        };

        // Sort command names for stable output
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

                // Variant header
                std::cout << "|   - Variant #" << (vi+1) << " — When: "
                          << (v.has_condition_ ? cond_to_string(v.condition_, "self", "parent", cond_to_string) : "(none)") << "\n";

                // Variant description block (placed BETWEEN header and Data-Layout)
                std::cout << "|     Description:\n";
                std::cout << "|       " << (v.doc_.empty() ? "-" : v.doc_) << "\n";

                // Data-Layout
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
};

} // namespace dsl
} // namespace fem
