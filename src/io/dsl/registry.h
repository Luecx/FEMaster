/**
 * @file registry.h
 * @brief Defines the command registry and DSL documentation helpers.
 *
 * `Registry` owns the complete syntax specification for one input dialect.
 * Commands are registered once and remain independent of semantic processing
 * order. `DeckParser` consumes the registry to validate and structure the source
 * deck, while FEMaster readers later execute selected parsed occurrences.
 *
 * The remaining methods render the registered grammar for command, variant,
 * token and search documentation. No parser-stage or command-activation state is
 * stored in the registry.
 *
 * @see Command
 * @see DeckParser
 *
 * @author Finn Eggers
 * @date 26.08.2026
 */

#pragma once

#include "../../core/logging.h"
#include "command.h"
#include "condition.h"
#include "pattern.h"
#include "pattern_element.h"
#include "segment.h"

#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace fem::io::dsl {

/**
 * @brief Owns and exposes all command specifications of one DSL dialect.
 */
struct Registry {
    std::unordered_map<std::string, Command> _map;

    /**
     * @brief Creates or retrieves a command and configures it in place.
     */
    void command(const std::string& name, const std::function<void(Command&)>& callback) {
        auto it = _map.find(name);
        if (it == _map.end()) it = _map.emplace(name, Command{name}).first;
        callback(it->second);
    }

    Command* find(const std::string& name) {
        const auto it = _map.find(name);
        return it == _map.end() ? nullptr : &it->second;
    }

    const Command* find(const std::string& name) const {
        const auto it = _map.find(name);
        return it == _map.end() ? nullptr : &it->second;
    }

    // ---------------------------------------------------------------------
    // Documentation helpers
    // ---------------------------------------------------------------------

    void print_index() const {
        std::vector<std::string> names;
        names.reserve(_map.size());
        for (const auto& entry : _map) names.push_back(entry.first);
        std::sort(names.begin(), names.end());

        std::cout << "= FEM DSL — Command Index =\n\n";
        for (const std::string& name : names) std::cout << name << "\n";
    }

    /**
     * @brief Prints the complete or compact grammar of one or all commands.
     */
    void print_help(const std::string& filter = {}, bool compact = false) const {
        const auto condition_to_string = [](const Condition& condition,
                                            const std::string& self,
                                            const std::string& parent,
                                            const auto& recurse) -> std::string {
            using Kind = Condition::Kind;

            switch (condition.kind) {
                case Kind::Always:
                    return "(none)";

                case Kind::ParentCommandIs: {
                    std::ostringstream stream;
                    stream << parent << ".command in { ";
                    for (std::size_t i = 0; i < condition.values.size(); ++i) {
                        if (i) stream << ", ";
                        stream << condition.values[i];
                    }
                    stream << " }";
                    return stream.str();
                }

                case Kind::ParentKeyEquals: {
                    std::ostringstream stream;
                    stream << parent << ".keys[\"" << condition.key << "\"] in { ";
                    for (std::size_t i = 0; i < condition.values.size(); ++i) {
                        if (i) stream << ", ";
                        stream << condition.values[i];
                    }
                    stream << " }";
                    return stream.str();
                }

                case Kind::ParentHasKey:
                    return parent + ".has_key(\"" + condition.key + "\")";

                case Kind::KeyPresent:
                    return self + ".has_key(\"" + condition.key + "\")";

                case Kind::KeyEquals: {
                    std::ostringstream stream;
                    stream << self << ".keys[\"" << condition.key << "\"] in { ";
                    for (std::size_t i = 0; i < condition.values.size(); ++i) {
                        if (i) stream << ", ";
                        stream << condition.values[i];
                    }
                    stream << " }";
                    return stream.str();
                }

                case Kind::ParentKeyBool:
                    return parent + ".keys[\"" + condition.key + "\"] == "
                         + (condition.bool_value ? "true" : "false");

                case Kind::KeyBool:
                    return self + ".keys[\"" + condition.key + "\"] == "
                         + (condition.bool_value ? "true" : "false");

                case Kind::All: {
                    std::ostringstream stream;
                    for (std::size_t i = 0; i < condition.children.size(); ++i) {
                        if (i) stream << " AND ";
                        stream << recurse(condition.children[i], self, parent, recurse);
                    }
                    return stream.str();
                }

                case Kind::Any: {
                    std::ostringstream stream;
                    for (std::size_t i = 0; i < condition.children.size(); ++i) {
                        if (i) stream << " OR ";
                        stream << recurse(condition.children[i], self, parent, recurse);
                    }
                    return stream.str();
                }

                case Kind::Not:
                    return "NOT(" + (condition.children.empty()
                        ? std::string("(none)")
                        : recurse(condition.children.front(), self, parent, recurse)) + ")";
            }

            return {};
        };

        const auto format_range = [](std::size_t minimum, std::size_t maximum) {
            std::ostringstream stream;
            stream << "[" << minimum << "..";
            if (maximum == std::numeric_limits<std::size_t>::max()) stream << "∞";
            else stream << maximum;
            stream << "]";
            return stream.str();
        };

        const auto render_pattern = [compact](const Pattern& pattern) {
            std::vector<std::string> headers;
            headers.reserve(pattern.required_tokens());
            std::size_t unnamed = 1;

            struct DocRow {
                std::string left;
                std::string doc;
                std::string type;
            };
            std::vector<DocRow> docs;

            for (const auto& element : pattern._elems) {
                const std::size_t count = element->count();
                const std::string base  = element->name_base();
                const std::string doc   = element->description();
                const char* type        = element->type_name();

                if (!base.empty()) {
                    if (count > 1) {
                        for (std::size_t i = 1; i <= count; ++i) {
                            headers.push_back(base + std::to_string(i));
                        }
                        docs.push_back({base + "1 - " + base + std::to_string(count), doc, type});
                    } else {
                        headers.push_back(base);
                        docs.push_back({base, doc, type});
                    }
                } else {
                    for (std::size_t i = 0; i < count; ++i) {
                        const std::string name = "t" + std::to_string(unnamed++);
                        headers.push_back(name);
                        docs.push_back({name, doc, type});
                    }
                }
            }

            if (!headers.empty()) {
                std::ostringstream stream;
                for (std::size_t i = 0; i < headers.size(); ++i) {
                    if (i) stream << ", ";
                    stream << headers[i];
                }
                fem::logging::info(true, "        ", stream.str());
            }

            if (compact) return;
            for (const DocRow& row : docs) {
                fem::logging::info(true,
                    "        ", std::left, std::setw(14), row.left + ":",
                    " ", row.doc.empty() ? "-" : row.doc, " [", row.type, "]");
            }
        };

        std::vector<std::string> names;
        names.reserve(_map.size());
        for (const auto& entry : _map) names.push_back(entry.first);
        std::sort(names.begin(), names.end());

        fem::logging::info(true, "= FEM DSL — Kurzreferenz =\n");

        for (const std::string& name : names) {
            if (!filter.empty() && name != filter) continue;
            const Command& command_spec = _map.at(name);

            fem::logging::info(true, "*", name);
            fem::logging::info(true, "Description:");
            fem::logging::info(true, "  ", command_spec.doc_.empty() ? std::string("-") : command_spec.doc_);

            if (command_spec.has_keyword_spec_) {
                const auto& entries = command_spec.keyword_spec_.entries();
                if (!entries.empty()) {
                    std::vector<const KeywordSpec::Entry*> keyword_entries;
                    keyword_entries.reserve(entries.size());
                    for (const auto& entry : entries) keyword_entries.push_back(&entry.second);
                    std::sort(keyword_entries.begin(), keyword_entries.end(),
                        [](const KeywordSpec::Entry* a, const KeywordSpec::Entry* b) {
                            return a->canonical < b->canonical;
                        });

                    fem::logging::info(true, "Keywords:");
                    for (const KeywordSpec::Entry* entry : keyword_entries) {
                        std::ostringstream header;
                        header << "  - " << entry->canonical;
                        if (entry->is_flag) header << " (flag)";
                        else header << " (" << (entry->required ? "required" : "optional") << ")";

                        if (entry->has_default) header << ", default=\"" << entry->default_value << "\"";
                        if (!entry->allowed.empty()) {
                            header << ", allowed={";
                            for (std::size_t i = 0; i < entry->allowed.size(); ++i) {
                                if (i) header << ", ";
                                header << entry->allowed[i];
                            }
                            header << "}";
                        }
                        fem::logging::info(true, header.str());
                        fem::logging::info(true, "    ", entry->doc.empty() ? std::string("-") : entry->doc);

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
            fem::logging::info(true, "  ",
                condition_to_string(command_spec.admit_, "self", "parent", condition_to_string));

            fem::logging::info(true, "Variants:");
            if (command_spec.variants_.empty()) {
                fem::logging::info(true, "  (none)");
                fem::logging::info(true);
                continue;
            }

            for (std::size_t variant_index = 0; variant_index < command_spec.variants_.size(); ++variant_index) {
                const Variant& variant = command_spec.variants_[variant_index];
                fem::logging::info(true,
                    "  - Variant #", variant_index + 1, " — When: ",
                    variant.has_condition_
                        ? condition_to_string(variant.condition_, "self", "parent", condition_to_string)
                        : std::string("(none)"));
                fem::logging::info(true, "    Description:");
                fem::logging::info(true, "      ", variant.doc_.empty() ? std::string("-") : variant.doc_);
                fem::logging::info(true, "    Data-Layout:");

                for (const Segment& segment : variant._segments) {
                    fem::logging::info(true,
                        "      • Lines ", format_range(segment._range.min_, segment._range.max_),
                        " (", segment._pattern.is_multiline() ? "multiline" : "single-line", "):");
                    render_pattern(segment._pattern);
                }
            }

            fem::logging::info(true);
        }
    }

    void print_variants(const std::string& command_name) const {
        const Command* command_spec = find(command_name);
        if (!command_spec) {
            std::cout << "Command not found: " << command_name << "\n";
            return;
        }

        const auto condition_to_string = [](const Condition& condition,
                                            const std::string& self,
                                            const std::string& parent,
                                            const auto& recurse) -> std::string {
            using Kind = Condition::Kind;
            switch (condition.kind) {
                case Kind::Always: return "(none)";
                case Kind::ParentCommandIs: {
                    std::ostringstream stream;
                    stream << parent << ".command in { ";
                    for (std::size_t i = 0; i < condition.values.size(); ++i) {
                        if (i) stream << ", ";
                        stream << condition.values[i];
                    }
                    stream << " }";
                    return stream.str();
                }
                case Kind::ParentKeyEquals: {
                    std::ostringstream stream;
                    stream << parent << ".keys[\"" << condition.key << "\"] in { ";
                    for (std::size_t i = 0; i < condition.values.size(); ++i) {
                        if (i) stream << ", ";
                        stream << condition.values[i];
                    }
                    stream << " }";
                    return stream.str();
                }
                case Kind::ParentHasKey: return parent + ".has_key(\"" + condition.key + "\")";
                case Kind::KeyPresent: return self + ".has_key(\"" + condition.key + "\")";
                case Kind::KeyEquals: {
                    std::ostringstream stream;
                    stream << self << ".keys[\"" << condition.key << "\"] in { ";
                    for (std::size_t i = 0; i < condition.values.size(); ++i) {
                        if (i) stream << ", ";
                        stream << condition.values[i];
                    }
                    stream << " }";
                    return stream.str();
                }
                case Kind::ParentKeyBool:
                    return parent + ".keys[\"" + condition.key + "\"] == "
                         + (condition.bool_value ? "true" : "false");
                case Kind::KeyBool:
                    return self + ".keys[\"" + condition.key + "\"] == "
                         + (condition.bool_value ? "true" : "false");
                case Kind::All: {
                    std::ostringstream stream;
                    for (std::size_t i = 0; i < condition.children.size(); ++i) {
                        if (i) stream << " AND ";
                        stream << recurse(condition.children[i], self, parent, recurse);
                    }
                    return stream.str();
                }
                case Kind::Any: {
                    std::ostringstream stream;
                    for (std::size_t i = 0; i < condition.children.size(); ++i) {
                        if (i) stream << " OR ";
                        stream << recurse(condition.children[i], self, parent, recurse);
                    }
                    return stream.str();
                }
                case Kind::Not:
                    return "NOT(" + (condition.children.empty()
                        ? std::string("(none)")
                        : recurse(condition.children.front(), self, parent, recurse)) + ")";
            }
            return {};
        };

        std::cout << "= Variants — " << command_name << " =\n\n";
        if (command_spec->variants_.empty()) {
            std::cout << "(none)\n";
            return;
        }

        for (std::size_t index = 0; index < command_spec->variants_.size(); ++index) {
            const Variant& variant = command_spec->variants_[index];
            std::cout << "- Variant #" << index + 1 << " — When: "
                      << (variant.has_condition_
                          ? condition_to_string(variant.condition_, "self", "parent", condition_to_string)
                          : "(none)") << "\n";
            std::cout << "  Doc: " << (variant.doc_.empty() ? "-" : variant.doc_) << "\n";
        }
    }

    void print_tokens(const std::string& command_name) const {
        const Command* command_spec = find(command_name);
        if (!command_spec) {
            std::cout << "Command not found: " << command_name << "\n";
            return;
        }

        const auto render_pattern = [](const Pattern& pattern, std::ostream& stream) {
            std::vector<std::string> headers;
            headers.reserve(pattern.required_tokens());
            std::size_t unnamed = 1;

            struct DocRow {
                std::string left;
                std::string doc;
                std::string type;
            };
            std::vector<DocRow> docs;

            for (const auto& element : pattern._elems) {
                const std::size_t count = element->count();
                const std::string base  = element->name_base();
                const std::string doc   = element->description();
                const char* type        = element->type_name();

                if (!base.empty()) {
                    if (count > 1) {
                        for (std::size_t i = 1; i <= count; ++i) headers.push_back(base + std::to_string(i));
                        docs.push_back({base + "1 - " + base + std::to_string(count), doc, type});
                    } else {
                        headers.push_back(base);
                        docs.push_back({base, doc, type});
                    }
                } else {
                    for (std::size_t i = 0; i < count; ++i) {
                        const std::string name = "t" + std::to_string(unnamed++);
                        headers.push_back(name);
                        docs.push_back({name, doc, type});
                    }
                }
            }

            stream << "Header: ";
            for (std::size_t i = 0; i < headers.size(); ++i) {
                if (i) stream << ", ";
                stream << headers[i];
            }
            stream << "\n";

            for (const DocRow& row : docs) {
                stream << "  - " << std::left << std::setw(12) << row.left + ":"
                       << " " << (row.doc.empty() ? "-" : row.doc) << " [" << row.type << "]\n";
            }
        };

        std::cout << "= Tokens — " << command_name << " =\n\n";
        if (command_spec->variants_.empty()) {
            std::cout << "(none)\n";
            return;
        }

        for (std::size_t variant_index = 0; variant_index < command_spec->variants_.size(); ++variant_index) {
            const Variant& variant = command_spec->variants_[variant_index];
            std::cout << "Variant #" << variant_index + 1 << ":\n";
            for (const Segment& segment : variant._segments) render_pattern(segment._pattern, std::cout);
            std::cout << "\n";
        }
    }

    void print_search(const std::string& query, bool regex = false) const {
        std::cout << "= Search: " << query << " =\n\n";

        std::vector<std::string> names;
        names.reserve(_map.size());
        for (const auto& entry : _map) names.push_back(entry.first);
        std::sort(names.begin(), names.end());

        std::regex expression;
        if (regex) expression = std::regex(query, std::regex::icase);

        int hits = 0;
        for (const std::string& name : names) {
            const Command& command_spec = _map.at(name);
            const auto matches = [&](const std::string& text) {
                if (regex) return std::regex_search(text, expression);

                std::string haystack = text;
                std::string needle   = query;
                std::transform(haystack.begin(), haystack.end(), haystack.begin(), ::tolower);
                std::transform(needle.begin(), needle.end(), needle.begin(), ::tolower);
                return haystack.find(needle) != std::string::npos;
            };

            if (matches(name) || matches(command_spec.doc_)) {
                std::cout << name << "  —  " << (command_spec.doc_.empty() ? "-" : command_spec.doc_) << "\n";
                ++hits;
            }
        }

        if (hits == 0) std::cout << "(no matches)\n";
    }

    void print_where_token(const std::string& token_substr) const {
        std::cout << "= Where Token: " << token_substr << " =\n\n";

        std::vector<std::string> names;
        names.reserve(_map.size());
        for (const auto& entry : _map) names.push_back(entry.first);
        std::sort(names.begin(), names.end());

        const auto contains_ci = [](const std::string& text, const std::string& substring) {
            std::string value = text;
            std::string query = substring;
            std::transform(value.begin(), value.end(), value.begin(), ::tolower);
            std::transform(query.begin(), query.end(), query.begin(), ::tolower);
            return value.find(query) != std::string::npos;
        };

        int hits = 0;
        for (const std::string& name : names) {
            const Command& command_spec = _map.at(name);
            bool found = false;

            for (const Variant& variant : command_spec.variants_) {
                for (const Segment& segment : variant._segments) {
                    for (const auto& element : segment._pattern._elems) {
                        if (contains_ci(element->name_base(), token_substr)
                         || contains_ci(element->description(), token_substr)) {
                            found = true;
                            break;
                        }
                    }
                    if (found) break;
                }
                if (found) break;
            }

            if (found) {
                std::cout << name << "\n";
                ++hits;
            }
        }

        if (hits == 0) std::cout << "(no matches)\n";
    }
};

} // namespace fem::io::dsl
