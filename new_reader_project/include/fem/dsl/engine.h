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
 *          - Choose the first matching **variant** (condition true or no condition).
 *          - Execute its **segments** according to their `LineRange` and `Pattern`.
 *          - Push this command as the new parent (scope goes “down”).
 *
 * ## Segment data acquisition semantics
 *
 *  - **Non-multiline pattern** (`pattern.is_multiline() == false`):
 *      * The pattern is a **per-line record layout**. Each DATA line is parsed independently,
 *        and the bound callback is invoked **once per line**.
 *      * The number of processed records (lines) must satisfy `min_ <= count <= max_`.
 *      * Early termination at `KEYWORD_LINE`/`END_OF_FILE` is allowed only if at least `min_`
 *        records have been processed; otherwise it is an error.
 *      * Excess tokens on a line are truncated to `required_tokens()`.
 *      * Missing/empty tokens on a line are tolerated only if the pattern's last element(s)
 *        configured `on_missing(...)` / `on_empty(...)`. This is applied via
 *        `Pattern::normalize_and_complete_tokens(...)`.
 *
 *  - **Multiline pattern** (`pattern.is_multiline() == true`):
 *      * The pattern defines a **record** that may span **multiple DATA lines**. The engine
 *        **aggregates tokens across consecutive data lines as needed** to satisfy the pattern’s
 *        `required_tokens()`. As soon as the quota is met (or can be met by appending defaults),
 *        the record is complete and the binder is invoked **once** for that record.
 *      * The segment then immediately attempts to read the **next record** in the same manner,
 *        repeating until a `KEYWORD_LINE` or `END_OF_FILE` is reached.
 *      * `LineRange.max_` is interpreted as a **per-record hard cap** on how many lines may be
 *        consumed. If a boundary is hit (KEYWORD/EOF) with a **partial record**:
 *          - If the pattern allows filling the remainder via `on_missing(...)`, the record is
 *            completed **and accepted**.
 *          - Otherwise, this is an error.
 *      * Note: `LineRange.min_` is intentionally **not** used as a stopping condition for
 *        multiline; the "read as many lines as needed to fill the pattern" rule governs completion.
 *
 *  - Comments and empty lines inside segments are ignored and do not count toward line limits.
 *
 * ## Lookahead
 *  - A one-line lookahead buffer preserves a `KEYWORD_LINE` when a segment ends exactly at a
 *    keyword boundary, so the next command is not lost.
 *
 * ## Errors
 *  - DATA line outside of an active command → error.
 *  - New KEYWORD while a non-multiline segment has processed fewer than `min_` records → error.
 *  - New KEYWORD/EOF while a multiline record cannot be completed (and defaults are insufficient)
 *    → error.
 *  - Unadmitted command in current scope → error (with scope trace).
 *  - No matching variant or missing segment binder → error.
 *
 * @see registry.h
 * @see file.h
 * @see line.h
 * @see keys.h
 */

#pragma once
#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>

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
     *      c. Select the first matching variant (condition true or no condition).
     *      d. Execute its segments using the rules described in the file header.
     *      e. Push the processed command as the new parent.
     *
     * @param file Input file object providing normalized lines (with include support).
     *
     * @throws std::runtime_error on structural or parsing inconsistencies.
     */
    void run(File& file) {
        using LT = LineType;

        // Scope: a stack of parent contexts (top = direct parent)
        std::vector<ParentInfo> scope;
        scope.push_back( ParentInfo{ "ROOT", Keys{} } );

        // One-line lookahead to preserve upcoming KEYWORD lines between segments/commands.
        Line buffered_line;
        bool have_buffered_line = false;

        /**
         * @brief Lambda: fetches the next line, honoring lookahead and optionally skipping ignorable lines.
         */
        auto pull_line = [&](bool skip_ignorable) -> Line {
            if (have_buffered_line) {
                have_buffered_line = false;
                return buffered_line;
            }
            if (skip_ignorable) {
                Line ln = file.next_line();
                return ln;
            }
            Line ln = file.next();
            return ln;
        };

        while (true) {
            Line ln = pull_line(true);

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

            // Lookup command in the registry
            const Command* spec = _reg.find(cmd);

            // Not registered → treat as structural parent and push; continue
            if (!spec) {
                scope.push_back( ParentInfo{ cmd, self_keys } );
                continue;
            }

            // Registered: attempt admission; climb up the scope until allowed
            int chosen_parent_index = admit_command(*spec, scope, self_keys);

            if (chosen_parent_index < 0) {
                // No admissible ancestor found → error with scope trace
                throw_unadmitted(cmd, scope);
            }

            // Trim scope to the admitted ancestor
            scope.resize(static_cast<std::size_t>(chosen_parent_index + 1));

            // Select the first matching variant (condition true or no condition)
            const Variant* chosen = select_variant(*spec, scope.back(), self_keys);
            if (!chosen)
                throw std::runtime_error("No matching variant for command: " + cmd);

            // Execute segments in order
            for (const auto& seg : chosen->_segments) {
                if (!seg._invoke)
                    throw std::runtime_error("No binder for a segment of " + cmd);

                const bool multiline = seg._pattern.is_multiline();
                if (multiline) {
                    read_multiline_segment(cmd,
                                           seg,
                                           pull_line,
                                           buffered_line,
                                           have_buffered_line);
                } else {
                    read_singleline_segment(cmd,
                                            seg,
                                            pull_line,
                                            buffered_line,
                                            have_buffered_line);
                }
            }

            // After successful processing, this keyword becomes the new parent (scope descends)
            scope.push_back( ParentInfo{ cmd, self_keys } );
        }
    }

private:
    /**
     * @brief Throws a detailed "not admitted" error including the current scope trace.
     *
     * @param cmd Command name that failed admission.
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
     * @param spec Command specification to check.
     * @param scope Current scope stack (top = direct parent).
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
     * @brief Chooses the first variant whose condition holds (or has no condition).
     *
     * @param spec Command specification.
     * @param parent The currently admitted parent context.
     * @param self_keys Keys from the current command's keyword line.
     * @return Pointer to the chosen variant, or nullptr if none matched.
     */
    const Variant* select_variant(const Command& spec,
                                  const ParentInfo& parent,
                                  const Keys& self_keys) const {
        for (const auto& v : spec.variants_) {
            if (!v.has_condition_ || v.condition_.eval(parent, self_keys)) {
                return &v;
            }
        }
        return nullptr;
    }

    /**
     * @brief Executes a non-multiline segment: one record per DATA line.
     *
     * The number of processed records (lines) must satisfy `min_ <= count <= max_`.
     * Each line is normalized and completed against the pattern (to honor on_empty/on_missing).
     *
     * @tparam PullFn Callable that yields a `Line`, honoring lookahead.
     * @param cmd Command name (for diagnostics).
     * @param seg Segment to execute.
     * @param pull_line Line puller bound to the current file and lookahead state.
     * @param buffered_line Lookahead buffer (to be set when we encounter the next keyword).
     * @param have_buffered_line Flag indicating whether the lookahead is active.
     *
     * @throws std::runtime_error On premature boundary before `min_` records or normalization failure.
     */
    template<class PullFn>
    void read_singleline_segment(const std::string& cmd,
                                 const Segment& seg,
                                 PullFn& pull_line,
                                 Line& buffered_line,
                                 bool& have_buffered_line) const {
        using LT = LineType;

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
                    buffered_line = dl;
                    have_buffered_line = true;
                    break;
                }
                throw std::runtime_error("Encountered next command while below minimum lines for segment of " + cmd);
            }

            if (dl.type() == LT::DATA_LINE) {
                std::vector<std::string> tokens;
                dl.append_values(tokens);

                // Truncate excess tokens
                if (tokens.size() > need_tokens) {
                    tokens.resize(need_tokens);
                }

                // Normalize/complete (apply on_empty/on_missing if configured)
                std::string err;
                if (!seg._pattern.normalize_and_complete_tokens(tokens, err)) {
                    if (err.empty()) err = "Token normalization failed for single-line segment";
                    throw std::runtime_error(err + " of " + cmd);
                }

                seg._invoke(tokens);
                ++records;
                continue;
            }

            // COMMENT / EMPTY → ignore and keep pulling
        }
    }

    /**
     * @brief Executes a multiline segment: repeating records, each aggregating multiple lines as needed.
     *
     * For each record:
     *  - Aggregate tokens across consecutive DATA lines until the pattern's token quota is met.
     *  - Respect `max_` as a per-record **line** cap if strictly required (only enforced when
     *    continuing to read another line would exceed the cap without being able to complete).
     *  - If a boundary (KEYWORD/EOF) is hit with a partial record:
     *      * Attempt to complete the record via `on_missing(...)`. If normalization succeeds,
     *        accept the record; otherwise, error.
     *  - Stop when a boundary appears **between records** (i.e., before reading any tokens for
     *    the next record). Preserve the boundary in the lookahead buffer if it is a keyword.
     *
     * @tparam PullFn Callable that yields a `Line`, honoring lookahead.
     * @param cmd Command name (for diagnostics).
     * @param seg Segment to execute.
     * @param pull_line Line puller bound to the current file and lookahead state.
     * @param buffered_line Lookahead buffer (to be set when we encounter the next keyword).
     * @param have_buffered_line Flag indicating whether the lookahead is active.
     *
     * @throws std::runtime_error On unfinishable partial records or hard cap violations.
     */
    template<class PullFn>
    void read_multiline_segment(const std::string& cmd,
                                const Segment& seg,
                                PullFn& pull_line,
                                Line& buffered_line,
                                bool& have_buffered_line) const {
        using LT = LineType;

        const std::size_t max_lines_per_record = seg._range.max_;
        const std::size_t need_tokens          = seg._pattern.required_tokens();

        while (true) {
            std::vector<std::string> tokens;
            std::size_t lines_in_record = 0;

            // Try to fill one record
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
                        // End between records → finish the segment
                        return;
                    }
                    // Partial record at EOF → try to finish via on_missing
                    std::string err;
                    if (!seg._pattern.normalize_and_complete_tokens(tokens, err)) {
                        if (err.empty()) err = "Unexpected EOF while filling multiline record";
                        throw std::runtime_error(err + " for " + cmd);
                    }
                    // Completed by defaults; invoke and finish segment at EOF
                    seg._invoke(tokens);
                    return;
                }

                if (dl.type() == LT::KEYWORD_LINE) {
                    if (tokens.empty()) {
                        // Boundary between records → stash and finish the segment
                        buffered_line = dl;
                        have_buffered_line = true;
                        return;
                    }
                    // Partial record at keyword → try to finish via on_missing
                    std::string err;
                    if (!seg._pattern.normalize_and_complete_tokens(tokens, err)) {
                        if (err.empty()) err = "Encountered next command before fulfilling multiline record";
                        throw std::runtime_error(err + " for " + cmd);
                    }
                    // Completed by defaults; invoke and stash the keyword for the outer loop
                    seg._invoke(tokens);
                    buffered_line = dl;
                    have_buffered_line = true;
                    return;
                }

                if (dl.type() == LT::DATA_LINE) {
                    dl.append_values(tokens);
                    ++lines_in_record;
                    continue;
                }

                // COMMENT / EMPTY → ignore and keep reading
            }

            // We have at least need_tokens; truncate excess if any
            if (tokens.size() > need_tokens) {
                tokens.resize(need_tokens);
            }

            // Normalize/complete (will be a no-op for fully filled, non-empty tokens)
            std::string err;
            if (!seg._pattern.normalize_and_complete_tokens(tokens, err)) {
                if (err.empty()) err = "Token normalization failed for multiline record";
                throw std::runtime_error(err + " of " + cmd);
            }

            // Invoke one complete record and loop to attempt reading the next
            seg._invoke(tokens);
        }
    }
};

} // namespace dsl
} // namespace fem
