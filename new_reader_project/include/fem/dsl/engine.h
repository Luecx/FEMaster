// ===== FILE: ./include/fem/dsl/engine.h =====
/**
 * @file
 * @brief High-level execution engine that interprets an input deck using the DSL registry.
 *
 * The engine streams normalized lines from `dsl::File`, maintains a scope stack of
 * parent commands, and executes registered commands by selecting admissible variants
 * and invoking their segments.
 *
 * ## Responsibilities
 *  - Maintain a scope stack of `(command, keys)` describing the current parent context.
 *  - For each encountered keyword line:
 *      * If the command is **not** registered: treat it as a structural parent and push it.
 *      * If the command **is** registered:
 *          - Attempt to **admit** it against the current scope (climb upwards until allowed).
 *          - Choose a **fitting variant** (preferring non-multiline variants first).
 *          - Execute its **segments** according to their `LineRange` and `Pattern`.
 *          - Push this command as the new parent (scope goes “down”).
 *
 * ## Variant selection (trial & backtrack)
 *  - Variants whose conditions hold are tried in two passes:
 *      1) All variants whose segments are **all non-multiline** (strict per-line layouts).
 *      2) All remaining variants (any **multiline** segment present).
 *  - Each candidate variant is **attempted** by reading lines, strictly validating
 *    token counts/types, and applying `on_empty`/`on_missing`:
 *      * If it **succeeds**, its effects are committed.
 *      * If it **fails**, all read lines are **rewound** (replayed) and the next variant
 *        is tried. The most informative error is surfaced if none fit.
 *
 * ## Segment data semantics
 *  - **Non-multiline pattern** (`pattern.is_multiline() == false`):
 *      * The pattern is a **per-line record layout**. Each DATA line yields **exactly one**
 *        record and the bound callback is invoked **once per line**.
 *      * **Strict token rule:** a DATA line is **invalid** if it contains **more**
 *        tokens than `required_tokens()` (no truncation). Fewer tokens are accepted
 *        only if `on_missing(...)` can complete the tail.
 *
 *  - **Multiline pattern** (`pattern.is_multiline() == true`):
 *      * The pattern defines a record that may span **multiple DATA lines**. Tokens are
 *        aggregated until the quota is met (or can be met via `on_missing(...)`).
 *      * The segment repeats to read additional records until a `KEYWORD_LINE` or `EOF`.
 *      * If a boundary is hit with a **partial record**, the engine tries to complete it
 *        via `on_missing(...)`. If that fails, the variant fails.
 *
 *  - Comments and empty lines inside segments are ignored and do not count toward line limits.
 *
 * ## Lookahead & rewind
 *  - A small **replay buffer** is used to push back lines when a variant attempt fails.
 *  - A one-line **keyword lookahead** (`buffered_line`) is propagated from a successful
 *    attempt so the outer loop sees the next command.
 *
 * ## Errors
 *  - DATA line outside of an active command → error.
 *  - Unadmitted command in current scope → error (with scope trace).
 *  - No variant fits → detailed error mentioning the closest attempt’s message.
 *
 * @see registry.h
 * @see file.h
 * @see line.h
 * @see keys.h
 */

#pragma once
#include <string>
#include <vector>
#include <deque>
#include <stdexcept>
#include <sstream>
#include <functional>

#include "registry.h"
#include "file.h"
#include "line.h"
#include "keys.h"

namespace fem {
namespace dsl {

/**
 * @class Engine
 * @brief Orchestrates reading, scope handling, variant selection, and segment execution.
 */
struct Engine {
    /**
     * @brief Reference to the command registry used for lookups and admission rules.
     */
    const Registry& _reg;

    /**
     * @brief Constructs an engine bound to a registry.
     *
     * @param r Command/DSL registry defining commands, variants, and segments.
     */
    explicit Engine(const Registry& r) : _reg(r) {}

    /**
     * @brief Executes the input deck by streaming from a file and applying the registry.
     *
     * Algorithm outline:
     *  1. Initialize scope with a synthetic `ROOT` parent (empty keys).
     *  2. Read normalized lines until `END_OF_FILE`.
     *  3. For each keyword line:
     *      a. If command is unregistered → push as structural parent and continue.
     *      b. If registered → climb scope upwards until the command is admitted.
     *      c. Try fitting variants (non-multiline first, then multiline) with backtracking.
     *      d. Execute its segments and push the processed command as the new parent.
     *
     * @param file Input file object providing normalized lines (with include support).
     *
     * @throws std::runtime_error on structural or parsing inconsistencies.
     */
    void run(File& file) {
        using LT = LineType;

        std::vector<ParentInfo> scope;
        scope.push_back( ParentInfo{ "ROOT", Keys{} } );

        // Global one-line keyword lookahead between commands (propagated from a successful attempt)
        Line  outer_buffered_kw;
        bool  outer_have_kw = false;

        // Replay buffer to "unread" lines (used only when a variant attempt fails)
        std::deque<Line> replay;

        // Base puller that honors replay and the outer keyword lookahead.
        auto base_pull = [&](bool skip_ignorable) -> Line {
            if (!replay.empty()) {
                Line x = replay.front();
                replay.pop_front();
                return x;
            }
            if (outer_have_kw) {
                outer_have_kw = false;
                return outer_buffered_kw;
            }
            return skip_ignorable ? file.next_line() : file.next();
        };

        while (true) {
            Line ln = base_pull(true);

            if (ln.type() == LT::END_OF_FILE)
                break;

            if (ln.type() == LT::DATA_LINE) {
                throw std::runtime_error("Unexpected DATA line without an active command context");
            }

            if (ln.type() != LT::KEYWORD_LINE)
                continue;

            // Extract command and normalized keys from the keyword line
            std::string cmd = ln.command();
            Keys self_keys  = Keys::from_keyword_line(ln);

            // Lookup command
            const Command* spec = _reg.find(cmd);

            // Not registered → structural node
            if (!spec) {
                scope.push_back( ParentInfo{ cmd, self_keys } );
                continue;
            }

            // Admission: climb scope
            int chosen_parent_index = admit_command(*spec, scope, self_keys);
            if (chosen_parent_index < 0) {
                throw_unadmitted(cmd, scope);
            }
            scope.resize(static_cast<std::size_t>(chosen_parent_index + 1));

            // Collect candidate variants whose condition holds
            std::vector<const Variant*> non_multi;
            std::vector<const Variant*> multi;
            for (const auto& v : spec->variants_) {
                if (v.has_condition_ && !v.condition_.eval(scope.back(), self_keys))
                    continue;
                bool any_multi = false;
                for (const auto& s : v._segments) {
                    if (s._pattern.is_multiline()) { any_multi = true; break; }
                }
                (any_multi ? multi : non_multi).push_back(&v);
            }

            // Two-pass attempt: non-multiline first, then multiline
            std::string last_err;
            bool        matched = false;

            auto try_pass = [&](const std::vector<const Variant*>& pass_list) {
                for (const Variant* vptr : pass_list) {
                    // Record every line consumed during this attempt
                    std::vector<Line> consumed;

                    // Local keyword lookahead inside the attempt
                    Line attempt_kw;
                    bool have_attempt_kw = false;

                    // Recording puller for the attempt (wraps base_pull)
                    auto attempt_pull = [&](bool skip_ignorable) -> Line {
                        Line got = base_pull(skip_ignorable);
                        consumed.push_back(got);
                        return got;
                    };

                    try {
                        execute_variant(cmd, *vptr, attempt_pull, attempt_kw, have_attempt_kw);
                        // Success → commit attempt's keyword lookahead to outer
                        if (have_attempt_kw) {
                            outer_buffered_kw = attempt_kw;
                            outer_have_kw = true;
                        }
                        matched = true;
                        return; // stop trying variants
                    } catch (const std::exception& e) {
                        last_err = e.what();

                        // Rewind: push consumed lines back to the front in reverse order
                        for (std::size_t i = consumed.size(); i-- > 0; ) {
                            replay.push_front(consumed[i]);
                        }
                        // Also clear any attempt-local keyword lookahead (not committed)
                        have_attempt_kw = false;
                    }
                }
            };

            try_pass(non_multi);
            if (!matched) try_pass(multi);

            if (!matched) {
                std::ostringstream os;
                os << "No variant of '" << cmd << "' fits the upcoming data";
                if (!last_err.empty()) os << ". Closest attempt failed with: " << last_err;
                throw std::runtime_error(os.str());
            }

            // After successful processing, this keyword becomes the new parent (scope descends)
            scope.push_back( ParentInfo{ cmd, self_keys } );
        }
    }

private:
    /**
     * @brief Throws a detailed "not admitted" error including the current scope trace.
     *
     * @param cmd   Command name that failed admission.
     * @param scope Current scope stack for diagnostics.
     *
     * @throws std::runtime_error Always throws with a composed message.
     */
    [[noreturn]] void throw_unadmitted(const std::string& cmd,
                                       const std::vector<ParentInfo>& scope) const {
        std::ostringstream os;
        os << "Command '" << cmd << "' not admitted in current scope (stack: ";
        for (std::size_t i = 0; i < scope.size(); ++i) {
            if (i) os << " > ";
            os << scope[i].command;
        }
        os << ")";
        throw std::runtime_error(os.str());
    }

    /**
     * @brief Attempts to admit a command by climbing the scope upwards.
     *
     * @param spec      Command specification to check.
     * @param scope     Current scope stack (top = direct parent).
     * @param self_keys Keys from the current command's keyword line.
     * @return Index of the chosen ancestor in the scope that admits the command, or -1 if none.
     */
    int admit_command(const Command& spec,
                      const std::vector<ParentInfo>& scope,
                      const Keys& self_keys) const {
        for (int i = static_cast<int>(scope.size()) - 1; i >= 0; --i) {
            const ParentInfo& candidate = scope[static_cast<std::size_t>(i)];
            if (spec.admit_.eval(candidate, self_keys)) {
                return i;
            }
        }
        return -1;
    }

    /**
     * @brief Executes all segments of a chosen variant using a provided puller.
     *
     * This function enforces:
     *  - Strict single-line token equality (no extra tokens allowed).
     *  - Multiline record aggregation with completion via on_missing if needed.
     *
     * @tparam PullFn     Callable that yields a `Line`, honoring the outer replay buffer.
     * @param cmd         Command name (for diagnostics).
     * @param var         Variant to execute.
     * @param pull_line   Line puller bound to the current file/buffers.
     * @param kw_out      Output slot for a stashed upcoming keyword (lookahead).
     * @param have_kw_out Output flag: true if `kw_out` is valid.
     */
    template<class PullFn>
    void execute_variant(const std::string& cmd,
                         const Variant& var,
                         PullFn& pull_line,
                         Line& kw_out,
                         bool& have_kw_out) const {
        using LT = LineType;

        for (const auto& seg : var._segments) {
            if (!seg._invoke)
                throw std::runtime_error("No binder for a segment of " + cmd);

            const bool multiline = seg._pattern.is_multiline();

            if (!multiline) {
                const std::size_t min_records = seg._range.min_;
                const std::size_t max_records = seg._range.max_;
                const std::size_t need_tokens = seg._pattern.required_tokens();

                std::size_t records = 0;

                while (true) {
                    if (records >= max_records) break;

                    Line dl = pull_line(false);

                    if (dl.type() == LT::END_OF_FILE) {
                        if (records >= min_records) break;
                        throw std::runtime_error("Unexpected EOF while reading single-line segment for " + cmd);
                    }

                    if (dl.type() == LT::KEYWORD_LINE) {
                        if (records >= min_records) {
                            kw_out = dl;
                            have_kw_out = true;
                            break;
                        }
                        throw std::runtime_error("Encountered next command while below minimum lines for segment of " + cmd);
                    }

                    if (dl.type() == LT::DATA_LINE) {
                        std::vector<std::string> tokens;
                        dl.append_values(tokens);

                        // STRICT: extra tokens invalidate the record (no truncation)
                        if (tokens.size() > need_tokens) {
                            std::ostringstream os;
                            os << "Too many tokens for single-line segment of " << cmd
                               << " (" << tokens.size() << " > " << need_tokens << ")";
                            throw std::runtime_error(os.str());
                        }

                        // Apply on_empty/on_missing normalization & completion
                        std::string err;
                        if (!seg._pattern.normalize_and_complete_tokens(tokens, err)) {
                            if (err.empty()) err = "Token normalization failed for single-line segment";
                            throw std::runtime_error(err + " of " + cmd);
                        }

                        // After normalization, we must have exactly need_tokens
                        if (tokens.size() != need_tokens) {
                            std::ostringstream os;
                            os << "Incomplete tokens for single-line segment of " << cmd
                               << " (" << tokens.size() << " != " << need_tokens << ")";
                            throw std::runtime_error(os.str());
                        }

                        seg._invoke(tokens);
                        ++records;
                        continue;
                    }

                    // COMMENT / EMPTY → ignore and keep pulling
                }

            } else {
                const std::size_t max_lines_per_record = seg._range.max_;
                const std::size_t need_tokens          = seg._pattern.required_tokens();

                while (true) {
                    std::vector<std::string> tokens;
                    std::size_t lines_in_record = 0;

                    // fill one record
                    while (tokens.size() < need_tokens) {
                        if (lines_in_record >= max_lines_per_record) {
                            std::ostringstream os;
                            os << "Exceeded maximum lines (" << max_lines_per_record
                               << ") before fulfilling multiline pattern for " << cmd;
                            throw std::runtime_error(os.str());
                        }

                        Line dl = pull_line(false);

                        if (dl.type() == LT::END_OF_FILE) {
                            if (tokens.empty()) {
                                // end between records
                                return;
                            }
                            // partial record at EOF → try to complete
                            std::string err;
                            if (!seg._pattern.normalize_and_complete_tokens(tokens, err)) {
                                if (err.empty()) err = "Unexpected EOF while filling multiline record";
                                throw std::runtime_error(err + " for " + cmd);
                            }
                            if (tokens.size() != need_tokens) {
                                std::ostringstream os;
                                os << "Unexpected EOF with incomplete multiline record for " << cmd;
                                throw std::runtime_error(os.str());
                            }
                            seg._invoke(tokens);
                            return;
                        }

                        if (dl.type() == LT::KEYWORD_LINE) {
                            if (tokens.empty()) {
                                // boundary between records
                                kw_out = dl;
                                have_kw_out = true;
                                return;
                            }
                            // partial record at keyword → try to complete
                            std::string err;
                            if (!seg._pattern.normalize_and_complete_tokens(tokens, err)) {
                                if (err.empty()) err = "Encountered next command before fulfilling multiline record";
                                throw std::runtime_error(err + " for " + cmd);
                            }
                            if (tokens.size() != need_tokens) {
                                std::ostringstream os;
                                os << "Next command arrived with incomplete multiline record for " << cmd;
                                throw std::runtime_error(os.str());
                            }
                            seg._invoke(tokens);
                            kw_out = dl;
                            have_kw_out = true;
                            return;
                        }

                        if (dl.type() == LT::DATA_LINE) {
                            dl.append_values(tokens);
                            ++lines_in_record;
                            continue;
                        }

                        // COMMENT / EMPTY → ignore
                    }

                    // We have ≥ need_tokens; **extra tokens are invalid** in multiline too
                    if (tokens.size() > need_tokens) {
                        std::ostringstream os;
                        os << "Too many tokens for multiline record of " << cmd
                           << " (" << tokens.size() << " > " << need_tokens << ")";
                        throw std::runtime_error(os.str());
                    }

                    // Normalize is still applied to honor on_empty (empties within lines)
                    std::string err;
                    if (!seg._pattern.normalize_and_complete_tokens(tokens, err)) {
                        if (err.empty()) err = "Token normalization failed for multiline record";
                        throw std::runtime_error(err + " of " + cmd);
                    }

                    if (tokens.size() != need_tokens) {
                        std::ostringstream os;
                        os << "Incomplete tokens for multiline record of " << cmd
                           << " (" << tokens.size() << " != " << need_tokens << ")";
                        throw std::runtime_error(os.str());
                    }

                    seg._invoke(tokens);
                    // loop to try another record; a boundary will be detected at the top
                }
            }
        }
    }
};

} // namespace dsl
} // namespace fem
