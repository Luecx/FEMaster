/**
 * @file deck_parser.h
 * @brief Parses a FEMaster DSL file once into a reusable `Deck`.
 *
 * `DeckParser` performs the complete syntactic work formerly repeated by every
 * semantic reader stage: keyword normalization, scope admission, variant
 * selection, segment validation, multiline aggregation and pattern completion.
 * Successful records are stored as `ParsedInvocation`s instead of being executed
 * immediately. The resulting deck can then be processed repeatedly in an
 * explicit dependency order without reopening or re-tokenizing the source file.
 *
 * Variant trial and rewind semantics are unchanged from the streaming engine.
 * The syntactic scope stack contains direct pointers to stable `ParsedCommand`
 * objects and starts at the deck's real ROOT node. Closing commands remain stored
 * as children of the scope they terminate while the scope stack is updated
 * immediately.
 *
 * @see Deck
 * @see ParsedCommand
 * @see Registry
 * @see Variant
 * @see Segment
 *
 * @author Finn Eggers
 * @date 26.08.2026
 */

#pragma once

#include "deck.h"
#include "file.h"
#include "keys.h"
#include "line.h"
#include "registry.h"

#include <algorithm>
#include <deque>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fem::io::dsl {

/**
 * @brief Converts one complete input file into a validated command tree.
 *
 * The parser owns no model state and never executes registered semantic
 * callbacks. Command and segment pointers stored in the resulting `Deck` refer
 * back to the supplied registry, which must therefore remain alive until deck
 * processing has completed.
 */
class DeckParser {
public:
    explicit DeckParser(const Registry& registry) : registry_(registry) {}

    /**
     * @brief Parses the complete file into a reusable deck.
     *
     * Scope admission follows the same ancestor-climbing rule as the previous
     * streaming engine. For every keyword, candidate variants are tried by
     * descending rank with full rewind on failure. The successful variant's
     * normalized invocations are retained below one `ParsedCommand` node.
     *
     * @param file Input source including nested INCLUDE expansion.
     * @return Fully validated parsed deck.
     */
    Deck parse(File& file) const {
        using LT = LineType;

        Deck deck;

        // Parsed nodes have stable addresses because every child is separately
        // owned. The scope therefore points directly into the resulting tree.
        std::vector<ParsedCommand*> scope{deck.root_.get()};

        auto pop_scope_to = [&](std::size_t new_size) {
            while (scope.size() > new_size) scope.pop_back();
        };

        // Preserve the existing one-keyword lookahead and variant rewind behavior.
        Line outer_buffered_keyword;
        bool outer_has_keyword = false;
        std::deque<Line> replay;

        auto base_pull = [&](bool skip_ignorable) -> Line {
            if (!replay.empty()) {
                Line line = replay.front();
                replay.pop_front();
                return line;
            }
            if (outer_has_keyword) {
                outer_has_keyword = false;
                return outer_buffered_keyword;
            }
            return skip_ignorable ? file.next_line() : file.next();
        };

        while (true) {
            Line line = base_pull(true);

            if (line.type() == LT::END_OF_FILE) break;
            if (line.type() == LT::DATA_LINE) {
                throw_at(line.location(), "Unexpected DATA line without an active command context");
            }
            if (line.type() != LT::KEYWORD_LINE) continue;

            // Resolve and validate the command keyword exactly once.
            const std::string command_name = line.command();
            Keys              self_keys    = Keys::from_keyword_line(line);
            const Command*    command      = registry_.find(command_name);

            if (!command) {
                throw_at(line.location(), "Unknown keyword '*" + command_name + "'");
            }

            if (command->has_keyword_spec_) {
                try {
                    self_keys.apply_spec(command->keyword_spec_, command_name);
                } catch (const std::exception& e) {
                    throw_at(line.location(), e.what());
                }
            }

            // Climb to the nearest syntactic ancestor admitting this command.
            const int parent_index = admit_command(*command, scope, self_keys);
            if (parent_index < 0) throw_unadmitted(line, scope);
            pop_scope_to(static_cast<std::size_t>(parent_index + 1));

            ParsedCommand* parent = scope.back();
            const ParentInfo parent_info = command_info(*parent);

            // Try all condition-compatible variants by rank with full input rewind.
            struct Candidate {
                const Variant* variant = nullptr;
                std::size_t    order   = 0;
            };

            std::vector<Candidate> candidates;
            candidates.reserve(command->variants_.size());
            for (std::size_t index = 0; index < command->variants_.size(); ++index) {
                const Variant& variant = command->variants_[index];
                if (variant.has_condition_ && !variant.condition_.eval(parent_info, self_keys)) continue;
                candidates.push_back(Candidate{&variant, index});
            }
            std::stable_sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
                if (a.variant->rank_ != b.variant->rank_) return a.variant->rank_ > b.variant->rank_;
                return a.order < b.order;
            });

            std::string                   last_error;
            bool                          matched = false;
            std::vector<ParsedInvocation> invocations;

            for (const Candidate& candidate : candidates) {
                std::vector<ParsedInvocation> candidate_invocations;
                std::vector<Line>             consumed;
                Line                          candidate_keyword;
                bool                          candidate_has_keyword = false;

                auto candidate_pull = [&](bool skip_ignorable) -> Line {
                    Line pulled = base_pull(skip_ignorable);
                    consumed.push_back(pulled);
                    return pulled;
                };

                try {
                    parse_variant(line, *candidate.variant, candidate_pull,
                                  candidate_keyword, candidate_has_keyword,
                                  candidate_invocations);

                    if (candidate_has_keyword) {
                        outer_buffered_keyword = candidate_keyword;
                        outer_has_keyword      = true;
                    }

                    invocations = std::move(candidate_invocations);
                    matched = true;
                    break;
                } catch (const std::exception& e) {
                    last_error = e.what();
                    for (std::size_t index = consumed.size(); index-- > 0;) {
                        replay.push_front(consumed[index]);
                    }
                }
            }

            if (!matched) {
                std::ostringstream message;
                message << "No variant of '" << command_name << "' fits the upcoming data";
                if (!last_error.empty()) message << ". Closest attempt failed with: " << last_error;
                throw_at(line.location(), message.str());
            }

            // A terminator may close any admitted scope except the real ROOT.
            if (command->closes_parent_ && parent == deck.root_.get()) {
                throw_at(line.location(), "Command '" + command_name + "' cannot close the ROOT scope");
            }

            // Every keyword is owned directly by its syntactic parent. The raw
            // pointer remains stable after insertion and can become the next scope.
            auto parsed = std::make_unique<ParsedCommand>();
            parsed->command_     = command;
            parsed->keys_        = std::move(self_keys);
            parsed->parent_      = parent;
            parsed->invocations_ = std::move(invocations);
            parsed->location_    = line.location();

            ParsedCommand* node = parsed.get();
            parent->children_.push_back(std::move(parsed));

            // Terminators stay in the tree but immediately close their parent.
            if (command->closes_parent_) {
                pop_scope_to(static_cast<std::size_t>(parent_index));
                continue;
            }

            // Ordinary commands become the new syntactic parent for later admission.
            scope.push_back(node);
        }

        return deck;
    }

private:
    const Registry& registry_;

    [[noreturn]] static void throw_at(const SourceLocation& location, const std::string& message) {
        const std::string source = location.str();
        throw std::runtime_error(source.empty() ? message : source + ": " + message);
    }

    static ParentInfo command_info(const ParsedCommand& command) {
        return ParentInfo{command.command_->name_, command.keys_};
    }

    [[noreturn]] static void throw_unadmitted(const Line& line,
                                              const std::vector<ParsedCommand*>& scope) {
        std::ostringstream message;
        message << "Command '" << line.command() << "' not admitted in current scope (stack: ";
        for (std::size_t index = 0; index < scope.size(); ++index) {
            if (index) message << " > ";
            message << scope[index]->command_->name_;
        }
        message << ")";
        throw_at(line.location(), message.str());
    }

    static int admit_command(const Command& command,
                             const std::vector<ParsedCommand*>& scope,
                             const Keys& self_keys) {
        for (int index = static_cast<int>(scope.size()) - 1; index >= 0; --index) {
            if (command.admit_.eval(command_info(*scope[static_cast<std::size_t>(index)]), self_keys)) {
                return index;
            }
        }
        return -1;
    }

    /**
     * @brief Parses every segment of one candidate variant without executing it.
     *
     * Single-line segments retain one invocation per DATA line. Multiline
     * patterns aggregate physical lines until one complete logical record has
     * been normalized. The implementation deliberately mirrors the previous
     * engine so existing `.inl` grammar definitions keep their exact behavior.
     */
    template<class PullFn>
    static void parse_variant(const Line& command_line,
                              const Variant& variant,
                              PullFn& pull_line,
                              Line& keyword_out,
                              bool& has_keyword_out,
                              std::vector<ParsedInvocation>& invocations) {
        using LT = LineType;

        const std::string& command_name = command_line.command();

        for (const Segment& segment : variant._segments) {
            const bool multiline = segment._pattern.is_multiline();

            if (!multiline) {
                const std::size_t min_records = segment._range.min_;
                const std::size_t max_records = segment._range.max_;
                const std::size_t need_tokens = segment._pattern.required_tokens();

                std::size_t records = 0;
                while (true) {
                    if (records >= max_records) break;

                    Line data_line = pull_line(false);

                    if (data_line.type() == LT::END_OF_FILE) {
                        if (records >= min_records) break;
                        throw_at(data_line.location(),
                            "Unexpected EOF while reading single-line segment for " + command_name);
                    }

                    if (data_line.type() == LT::KEYWORD_LINE) {
                        if (records >= min_records) {
                            keyword_out     = data_line;
                            has_keyword_out = true;
                            break;
                        }
                        throw_at(data_line.location(),
                            "Encountered next command while below minimum lines for segment of " + command_name);
                    }

                    if (data_line.type() == LT::DATA_LINE) {
                        std::vector<std::string> tokens;
                        data_line.append_values(tokens);

                        if (tokens.size() > need_tokens) {
                            std::ostringstream message;
                            message << "Too many tokens for single-line segment of " << command_name
                                    << " (" << tokens.size() << " > " << need_tokens << ")";
                            throw_at(data_line.location(), message.str());
                        }

                        std::string error;
                        if (!segment._pattern.normalize_and_complete_tokens(tokens, error)) {
                            if (error.empty()) error = "Token normalization failed for single-line segment";
                            throw_at(data_line.location(), error + " of " + command_name);
                        }

                        if (tokens.size() != need_tokens) {
                            std::ostringstream message;
                            message << "Incomplete tokens for single-line segment of " << command_name
                                    << " (" << tokens.size() << " != " << need_tokens << ")";
                            throw_at(data_line.location(), message.str());
                        }

                        if (segment._invoke) {
                            invocations.push_back(ParsedInvocation{&segment, std::move(tokens), data_line.location()});
                        }
                        ++records;
                    }
                }
                continue;
            }

            const std::size_t max_lines_per_record = segment._range.max_;
            const std::size_t need_tokens          = segment._pattern.required_tokens();

            while (true) {
                std::vector<std::string> tokens;
                std::size_t              lines_in_record = 0;
                SourceLocation           record_location;

                while (tokens.size() < need_tokens) {
                    if (lines_in_record >= max_lines_per_record) {
                        std::ostringstream message;
                        message << "Exceeded maximum lines (" << max_lines_per_record
                                << ") before fulfilling multiline pattern for " << command_name;
                        throw_at(record_location.source ? record_location : command_line.location(), message.str());
                    }

                    Line data_line = pull_line(false);

                    if (data_line.type() == LT::END_OF_FILE) {
                        if (tokens.empty()) return;

                        std::string error;
                        if (!segment._pattern.normalize_and_complete_tokens(tokens, error)) {
                            if (error.empty()) error = "Unexpected EOF while filling multiline record";
                            throw_at(record_location, error + " for " + command_name);
                        }
                        if (tokens.size() != need_tokens) {
                            throw_at(record_location,
                                "Unexpected EOF with incomplete multiline record for " + command_name);
                        }
                        if (segment._invoke) {
                            invocations.push_back(ParsedInvocation{&segment, std::move(tokens), record_location});
                        }
                        return;
                    }

                    if (data_line.type() == LT::KEYWORD_LINE) {
                        if (tokens.empty()) {
                            keyword_out     = data_line;
                            has_keyword_out = true;
                            return;
                        }

                        std::string error;
                        if (!segment._pattern.normalize_and_complete_tokens(tokens, error)) {
                            if (error.empty()) error = "Encountered next command before fulfilling multiline record";
                            throw_at(record_location, error + " for " + command_name);
                        }
                        if (tokens.size() != need_tokens) {
                            throw_at(record_location,
                                "Next command arrived with incomplete multiline record for " + command_name);
                        }
                        if (segment._invoke) {
                            invocations.push_back(ParsedInvocation{&segment, std::move(tokens), record_location});
                        }
                        keyword_out     = data_line;
                        has_keyword_out = true;
                        return;
                    }

                    if (data_line.type() == LT::DATA_LINE) {
                        if (lines_in_record == 0) record_location = data_line.location();
                        data_line.append_values(tokens);
                        ++lines_in_record;
                    }
                }

                if (tokens.size() > need_tokens) {
                    std::ostringstream message;
                    message << "Too many tokens for multiline record of " << command_name
                            << " (" << tokens.size() << " > " << need_tokens << ")";
                    throw_at(record_location, message.str());
                }

                std::string error;
                if (!segment._pattern.normalize_and_complete_tokens(tokens, error)) {
                    if (error.empty()) error = "Token normalization failed for multiline record";
                    throw_at(record_location, error + " of " + command_name);
                }
                if (tokens.size() != need_tokens) {
                    std::ostringstream message;
                    message << "Incomplete tokens for multiline record of " << command_name
                            << " (" << tokens.size() << " != " << need_tokens << ")";
                    throw_at(record_location, message.str());
                }

                if (segment._invoke) {
                    invocations.push_back(ParsedInvocation{&segment, std::move(tokens), record_location});
                }
            }
        }
    }
};

} // namespace fem::io::dsl
