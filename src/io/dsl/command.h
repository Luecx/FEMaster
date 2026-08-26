/**
 * @file command.h
 * @brief Declares `Command`, the top-level specification for a DSL keyword.
 *
 * A `Command` describes one registered keyword independently of when its
 * semantic callback will execute. It defines scope admission, keyword arguments,
 * candidate data variants and optional entry/exit hooks. `DeckParser` applies
 * the syntactic rules once and stores the resulting command occurrences;
 * higher-level readers later choose the semantic execution order explicitly.
 *
 * @see variant.h
 * @see condition.h
 * @see keyword.h
 * @see deck_parser.h
 *
 * @author Finn Eggers
 * @date 26.08.2026
 */

#pragma once

#include "condition.h"
#include "keyword.h"
#include "variant.h"

#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace fem::io::dsl {

/**
 * @brief Complete grammar and callback specification for one DSL keyword.
 *
 * The command owns no parser-stage state. Admission and variants describe only
 * syntax, while semantic timing remains a responsibility of the reader that
 * executes the parsed deck.
 */
struct Command {
    // Command identity, scope admission and data layouts
    std::string          name_;
    Condition            admit_ = Condition{};
    std::vector<Variant> variants_;
    std::string          doc_;

    // Scope behavior and semantic callbacks
    bool closes_parent_ = false;

    std::function<void(const ParentInfo&, const Keys&)> on_enter_;
    std::function<void(const Keys&)>                    on_exit_;

    // Keyword-argument normalization and validation
    KeywordSpec keyword_spec_;
    bool        has_keyword_spec_ = false;

    explicit Command(std::string name) : name_(std::move(name)) {}

    /**
     * @brief Restricts this keyword to parents satisfying the supplied condition.
     */
    Command& allow_if(Condition condition) {
        admit_ = std::move(condition);
        return *this;
    }

    /**
     * @brief Appends one candidate data layout for this command.
     *
     * `DeckParser` evaluates compatible variants in descending rank order and
     * retains the first layout whose complete upcoming data block fits.
     */
    Command& variant(Variant variant) {
        variants_.push_back(std::move(variant));
        return *this;
    }

    /**
     * @brief Sets the short description shown by DSL documentation output.
     */
    Command& doc(std::string description) {
        doc_ = std::move(description);
        return *this;
    }

    /**
     * @brief Marks this keyword as a terminator of its admitted parent scope.
     *
     * During parsing the terminator remains stored as a child occurrence of the
     * parent it closes, while the syntactic scope stack immediately returns to
     * the parent's parent. Its own callbacks can therefore still be executed
     * later at the explicit semantic point chosen by the reader.
     */
    Command& closes_parent() {
        closes_parent_ = true;
        return *this;
    }

    /**
     * @brief Registers a parent-aware hook executed before stored data records.
     */
    Command& on_enter(std::function<void(const ParentInfo&, const Keys&)> callback) {
        on_enter_ = std::move(callback);
        return *this;
    }

    /**
     * @brief Registers a key-only hook executed before stored data records.
     */
    Command& on_enter(std::function<void(const Keys&)> callback) {
        on_enter_ = [callback = std::move(callback)](const ParentInfo&, const Keys& keys) {
            callback(keys);
        };
        return *this;
    }

    /**
     * @brief Registers a hook executed when semantic processing leaves this scope.
     *
     * The reader decides when a parsed scope is left. This preserves command-local
     * finalization such as restoring the default Part or executing a completed
     * load case without coupling that behavior to source-file streaming.
     */
    Command& on_exit(std::function<void(const Keys&)> callback) {
        on_exit_ = std::move(callback);
        return *this;
    }

    /**
     * @brief Declares keyword arguments, aliases, defaults and value constraints.
     */
    Command& keyword(KeywordSpec spec) {
        keyword_spec_     = std::move(spec);
        has_keyword_spec_ = true;
        return *this;
    }
};

} // namespace fem::io::dsl
